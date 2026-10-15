proxlq_addon_root <- function(parent_root = project_root()) {
  normalizePath(file.path(parent_root, "proxlq_addon"), mustWork = TRUE)
}

source_proxlq_addon <- function(parent_root = project_root(), compile_cpp = TRUE) {
  parent_root <- normalizePath(parent_root, mustWork = TRUE)
  source(file.path(parent_root, "config", "simulation_config.R"), local = .GlobalEnv)
  for (file in c("dgp.R", "bridges.R", "bandwidth.R", "inference.R"))
    source(file.path(parent_root, "R", file), local = .GlobalEnv)
  if (isTRUE(compile_cpp)) {
    load_fast_cv(parent_root)
    addon_root <- proxlq_addon_root(parent_root)
    cache <- Sys.getenv("PDR_PROXLQ_CPP_CACHE", unset = "")
    if (!nzchar(cache)) cache <- file.path(addon_root, ".cache", "Rcpp")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    if (!exists("epan_derivative_lq_loocv_cpp", mode = "function", inherits = TRUE)) {
      if (!requireNamespace("Rcpp", quietly = TRUE))
        stop("Package 'Rcpp' is required for ProxLQ CV", call. = FALSE)
      Rcpp::sourceCpp(file.path(addon_root, "src", "fast_lq_cv.cpp"),
                      cacheDir = cache, rebuild = FALSE, showOutput = FALSE)
    }
  }
  invisible(parent_root)
}

proxlq_addon_config <- function(profile = c("production", "pilot", "test"),
                                parent_root = project_root()) {
  profile <- match.arg(profile)
  base <- simulation_config(profile)
  production_parent <- simulation_config("production")
  parent_id <- config_fingerprint(production_parent, parent_root)
  if (profile == "production")
    assert_true(grepl("^[0-9a-f]{8}$", parent_id),
                "Invalid packaged parent production fingerprint")
  list(
    version = "proxlq-derivative-addon-v1.0.0",
    profile = profile,
    parent_config_id = parent_id,
    dgp = base$dgp,
    sample_sizes = base$sample_sizes,
    repetitions = base$repetitions,
    base_seed = base$base_seed,
    folds = base$folds,
    curve_grid = base$curve_grid,
    point_grid = base$point_grid,
    alpha = base$alpha,
    scenarios = base$scenarios,
    methods = c("ProxLQ-US-log10", "ProxLQ-US-log2over7"),
    polynomial_degree = 2L,
    derivative_order = 1L,
    raw_selector = "CV1-LQ",
    cv1_multipliers = base$cv1_multipliers,
    cv1_refine = base$cv1_refine,
    derivative_rate_exponent = -1 / 7,
    derivative_cv_domain = base$derivative_cv_domain,
    derivative_cv_flat_domain = base$derivative_cv_flat_domain,
    base_cfg = base
  )
}

addon_config_fingerprint <- function(cfg, parent_root = project_root(),
                                     addon_root = proxlq_addon_root(parent_root)) {
  files <- sort(c(
    file.path(parent_root, "config", "simulation_config.R"),
    file.path(parent_root, "R", c("utils.R", "dgp.R", "bridges.R", "bandwidth.R",
                                   "inference.R")),
    file.path(parent_root, "src", "fast_cv.cpp"),
    list.files(file.path(addon_root, "R"), pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(addon_root, "src"), pattern = "[.]cpp$", full.names = TRUE),
    list.files(file.path(addon_root, "scripts"), pattern = "^(run_worker|validate_math)[.]R$",
               full.names = TRUE)
  ))
  manifest <- paste(substring(files, nchar(parent_root) + 2L),
                    unname(tools::md5sum(files)), collapse = "\n")
  payload <- paste(
    paste(capture.output(dput(cfg, control = c("keepNA", "keepInteger"))), collapse = "\n"),
    manifest, sep = "\n--SOURCE--\n"
  )
  raw <- as.integer(charToRaw(enc2utf8(payload)))
  hash <- 0
  for (value in raw) hash <- (hash * 131 + value) %% 2147483647
  sprintf("%08x", as.integer(hash))
}

validate_proxlq_design <- function(cfg, n = NULL) {
  assert_true(identical(cfg$methods,
                        c("ProxLQ-US-log10", "ProxLQ-US-log2over7")),
              "ProxLQ method list mismatch")
  assert_true(cfg$polynomial_degree == 2L && cfg$derivative_order == 1L,
              "ProxLQ polynomial specification mismatch")
  assert_true(identical(cfg$scenarios$scenario, c("CC", "CM", "MC", "MM")),
              "Scenario list mismatch")
  if (cfg$profile == "production")
    assert_true(length(cfg$curve_grid) == 61L && length(cfg$point_grid) == 9L,
                "Production evaluation-grid shape mismatch")
  assert_true(max(abs(range(cfg$curve_grid) - c(-0.9, 0.9))) < 1e-14,
              "Curve grid mismatch")
  assert_true(all(grid_membership(cfg$point_grid, cfg$curve_grid)),
              "Point grid is not a curve-grid subset")
  assert_true(cfg$derivative_rate_exponent == -1 / 7,
              "Derivative reference rate must be n^(-1/7)")
  if (!is.null(n))
    assert_true(n >= 5L * cfg$base_cfg$cv_min_neighbors,
                "n is too small for configured local fits")
  invisible(TRUE)
}

derivative_lq_loocv <- function(a, xi, inverse_density, bandwidths, cfg) {
  epan_derivative_lq_loocv_cpp(
    as.numeric(a), as.matrix(xi), as.numeric(inverse_density),
    sort(unique(as.numeric(bandwidths))),
    cfg$derivative_cv_domain[[1L]], cfg$derivative_cv_domain[[2L]],
    cfg$derivative_cv_flat_domain[[1L]], cfg$derivative_cv_flat_domain[[2L]],
    cfg$base_cfg$cv_min_neighbors
  )
}

select_proxlq_bandwidth <- function(a, xi, cfg) {
  density <- derivative_density(a, cfg$base_cfg)
  base <- reference_bandwidth(a, cfg$derivative_rate_exponent)
  initial <- base * cfg$cv1_multipliers
  first <- derivative_lq_loocv(a, xi, density$inverse, initial, cfg)
  allowed <- if (isTRUE(cfg$cv1_refine))
    adjacent_refinement(first$bandwidth, first$risk_lq) else
      rep(list(initial), ncol(xi))
  candidates <- sort(unique(unlist(allowed)))
  final <- derivative_lq_loocv(a, xi, density$inverse, candidates, cfg)
  selection <- select_from_allowed(final$bandwidth, final$risk_lq, allowed)
  list(base = base, initial = initial, density = density, final = final,
       selection = selection)
}

proxlq_divisor <- function(n, method) {
  if (method == "ProxLQ-US-log10") return(log10(n))
  if (method == "ProxLQ-US-log2over7") return(log(n)^(2 / 7))
  stop("Unknown ProxLQ undersmoothing method", call. = FALSE)
}

proxlq_candidate_allowed <- function(h, allowed) {
  vapply(h, function(value) any(abs(log(value / allowed)) < 1e-12), logical(1L))
}

evaluate_proxlq_pointwise <- function(data, crossfit, scenario, method, h, cfg) {
  scenario_row <- cfg$scenarios[cfg$scenarios$scenario == scenario, , drop = FALSE]
  assert_true(nrow(scenario_row) == 1L, "Unknown nuisance scenario")
  outcome_model <- scenario_row$outcome_model[[1L]]
  result <- eif_at_grid(
    data$A, crossfit$xi[, scenario], crossfit$fold,
    crossfit$h_delta[, outcome_model], cfg$curve_grid, h,
    cfg$polynomial_degree, cfg$derivative_order, cfg$base_cfg
  )
  z <- qnorm(1 - cfg$alpha / 2)
  curve <- data.frame(
    estimand = "derivative", scenario = scenario, method = method,
    a = cfg$curve_grid,
    on_point_grid = grid_membership(cfg$curve_grid, cfg$point_grid),
    truth = theta_derivative_truth(cfg$curve_grid),
    estimate = result$estimate, se = result$se,
    ci_lower = result$estimate - z * result$se,
    ci_upper = result$estimate + z * result$se,
    bandwidth = h, raw_bandwidth = NA_real_,
    polynomial_degree = cfg$polynomial_degree,
    derivative_order = cfg$derivative_order,
    phi_mean = result$phi_mean,
    max_condition = result$condition,
    reproduction_error = result$reproduction_error,
    local_n = result$local_n,
    stringsAsFactors = FALSE
  )
  curve$pointwise_coverage <- curve$truth >= curve$ci_lower &
    curve$truth <= curve$ci_upper
  curve$ci_width <- curve$ci_upper - curve$ci_lower
  curve$squared_error <- (curve$estimate - curve$truth)^2
  curve
}

proxlq_cv_diagnostics <- function(n, rep, cfg, bandwidth) {
  rows <- vector("list", nrow(cfg$scenarios))
  for (index in seq_len(nrow(cfg$scenarios))) {
    scenario <- cfg$scenarios$scenario[[index]]
    h <- bandwidth$final$bandwidth
    risk <- bandwidth$final$risk_lq[, scenario]
    selection <- bandwidth$selection
    allowed <- proxlq_candidate_allowed(h, selection$allowed[[scenario]])
    rows[[index]] <- data.frame(
      n = n, rep = rep, scenario = scenario, selector = cfg$raw_selector,
      bandwidth = h[allowed], risk = risk[allowed],
      finite_candidate = is.finite(risk[allowed]),
      selected = abs(log(h[allowed] / selection$h[[scenario]])) < 1e-12,
      invalid_evaluations = bandwidth$final$invalid_lq[allowed],
      target_evaluations = bandwidth$final$tapered_evaluations[allowed],
      min_neighbors = bandwidth$final$min_neighbors[allowed],
      min_pivot_ratio = bandwidth$final$min_pivot_ratio[allowed],
      density_floor_fraction = bandwidth$density$floor_fraction,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_one_proxlq_replication <- function(n, rep, cfg, parent_root = project_root()) {
  validate_proxlq_design(cfg, n)
  started <- proc.time()[[3L]]
  dgp_seed <- seed_for(n, rep, 1L, cfg$base_seed)
  fold_seed <- seed_for(n, rep, 2L, cfg$base_seed)
  data <- simulate_proximal_data(n, dgp_seed)
  fold <- balanced_folds(n, cfg$folds, fold_seed)

  nuisance_started <- proc.time()[[3L]]
  crossfit <- crossfit_bridges(data, fold, cfg$base_cfg)
  nuisance_seconds <- proc.time()[[3L]] - nuisance_started

  bandwidth_started <- proc.time()[[3L]]
  bandwidth <- select_proxlq_bandwidth(data$A, crossfit$xi, cfg)
  bandwidth_seconds <- proc.time()[[3L]] - bandwidth_started

  inference_started <- proc.time()[[3L]]
  curves <- vector("list", length(cfg$methods) * nrow(cfg$scenarios))
  bandwidth_rows <- vector("list", length(curves))
  row_index <- 0L
  scale <- robust_exposure_scale(data$A)
  for (scenario in cfg$scenarios$scenario) {
    raw_h <- bandwidth$selection$h[[scenario]]
    for (method in cfg$methods) {
      row_index <- row_index + 1L
      divisor <- proxlq_divisor(n, method)
      final_h <- raw_h / divisor
      curve <- evaluate_proxlq_pointwise(data, crossfit, scenario, method, final_h, cfg)
      curve$raw_bandwidth <- raw_h
      curves[[row_index]] <- curve
      bandwidth_rows[[row_index]] <- data.frame(
        n = n, rep = rep, scenario = scenario, method = method,
        raw_h = raw_h, divisor = divisor, h = final_h,
        normalized_raw_h = raw_h / (scale * n^(-1 / 7)),
        normalized_h = final_h / (scale * n^(-1 / 7)),
        n_h3 = n * final_h^3, n_h7 = n * final_h^7,
        selected_risk = bandwidth$selection$selected_risk[[scenario]],
        second_risk = bandwidth$selection$second_risk[[scenario]],
        risk_gap = bandwidth$selection$risk_gap[[scenario]],
        normalized_risk_gap = bandwidth$selection$normalized_risk_gap[[scenario]],
        lower_hit = bandwidth$selection$lower_hit[[scenario]],
        upper_hit = bandwidth$selection$upper_hit[[scenario]],
        boundary_hit = bandwidth$selection$boundary_hit[[scenario]],
        stringsAsFactors = FALSE
      )
    }
  }
  inference_seconds <- proc.time()[[3L]] - inference_started
  curve <- do.call(rbind, curves)
  curve$n <- n; curve$rep <- rep
  curve <- curve[c("n", "rep", setdiff(names(curve), c("n", "rep")))]
  rownames(curve) <- NULL

  list(
    metadata = list(
      version = cfg$version, profile = cfg$profile,
      config_id = addon_config_fingerprint(cfg, parent_root),
      parent_config_id = cfg$parent_config_id,
      n = n, rep = rep, dgp_seed = dgp_seed, fold_seed = fold_seed,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      session = capture.output(sessionInfo())
    ),
    curve = curve,
    bandwidth = do.call(rbind, bandwidth_rows),
    cv_diagnostics = proxlq_cv_diagnostics(n, rep, cfg, bandwidth),
    nuisance_diagnostics = crossfit$diagnostics,
    dgp_diagnostics = dgp_diagnostics(data, range(cfg$curve_grid)),
    timing = c(nuisance_seconds = nuisance_seconds,
               bandwidth_seconds = bandwidth_seconds,
               inference_seconds = inference_seconds,
               total_seconds = proc.time()[[3L]] - started)
  )
}

proxlq_result_path <- function(outdir, config_id, n, rep) {
  file.path(outdir, paste0("config_", config_id), paste0("n_", n),
            sprintf("rep_%04d.rds", rep))
}

validate_proxlq_part <- function(path, cfg, n, rep,
                                 parent_root = project_root()) {
  if (!file.exists(path)) return(FALSE)
  answer <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(answer) || is.null(answer$metadata)) return(FALSE)
  metadata_ok <- identical(answer$metadata$config_id,
                           addon_config_fingerprint(cfg, parent_root)) &&
    identical(answer$metadata$parent_config_id, cfg$parent_config_id) &&
    identical(as.integer(answer$metadata$n), as.integer(n)) &&
    identical(as.integer(answer$metadata$rep), as.integer(rep)) &&
    identical(as.integer(answer$metadata$dgp_seed),
              seed_for(n, rep, 1L, cfg$base_seed)) &&
    identical(as.integer(answer$metadata$fold_seed),
              seed_for(n, rep, 2L, cfg$base_seed))
  if (!metadata_ok || !is.data.frame(answer$curve) ||
      !is.data.frame(answer$bandwidth) || !is.data.frame(answer$cv_diagnostics))
    return(FALSE)

  curve <- answer$curve
  expected_rows <- nrow(cfg$scenarios) * length(cfg$methods) * length(cfg$curve_grid)
  key <- paste(curve$scenario, curve$method, sprintf("%.12g", curve$a), sep = "|")
  finite_columns <- c("truth", "estimate", "se", "ci_lower", "ci_upper",
                      "bandwidth", "raw_bandwidth", "phi_mean", "max_condition",
                      "reproduction_error", "local_n", "ci_width", "squared_error")
  curve_ok <- nrow(curve) == expected_rows && !anyDuplicated(key) &&
    all(curve$n == n) && all(curve$rep == rep) &&
    all(curve$estimand == "derivative") &&
    all(curve$polynomial_degree == 2L) && all(curve$derivative_order == 1L) &&
    all(is.finite(as.matrix(curve[finite_columns]))) && all(curve$se > 0) &&
    all(curve$local_n >= cfg$base_cfg$inference_min_neighbors) &&
    max(curve$reproduction_error) <= cfg$base_cfg$local_reproduction_tolerance

  bandwidth <- answer$bandwidth
  bandwidth_key <- paste(bandwidth$scenario, bandwidth$method, sep = "|")
  formula_ok <- vapply(seq_len(nrow(bandwidth)), function(i) {
    expected <- bandwidth$raw_h[[i]] / proxlq_divisor(n, bandwidth$method[[i]])
    abs(log(bandwidth$h[[i]] / expected)) < 1e-13
  }, logical(1L))
  bandwidth_ok <- nrow(bandwidth) == nrow(cfg$scenarios) * length(cfg$methods) &&
    !anyDuplicated(bandwidth_key) && all(formula_ok) &&
    all(is.finite(as.matrix(bandwidth[c("raw_h", "divisor", "h",
                                         "normalized_raw_h", "normalized_h",
                                         "n_h3", "n_h7")]))) && all(bandwidth$h > 0)

  cv <- answer$cv_diagnostics
  cv_group <- split(cv, cv$scenario)
  cv_ok <- length(cv_group) == nrow(cfg$scenarios) &&
    all(vapply(cv_group, function(value)
      sum(value$selected) == 1L && all(value$finite_candidate[value$selected]) &&
        all(is.finite(value$bandwidth)) && all(value$bandwidth > 0), logical(1L)))
  timing_ok <- length(answer$timing) == 4L && all(is.finite(answer$timing)) &&
    all(answer$timing >= 0)
  isTRUE(curve_ok && bandwidth_ok && cv_ok && timing_ok)
}

proxlq_failure_path <- function(path) sub("[.]rds$", ".failed.rds", path)

execute_proxlq_replication <- function(n, rep, cfg, parent_root, outdir,
                                       force = FALSE, worker = "single") {
  config_id <- addon_config_fingerprint(cfg, parent_root)
  path <- proxlq_result_path(outdir, config_id, n, rep)
  if (!force && validate_proxlq_part(path, cfg, n, rep, parent_root)) {
    return(data.frame(n = n, rep = rep, worker = worker, status = "skipped_valid",
                      seconds = 0, path = path, detail = "", stringsAsFactors = FALSE))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(path, ".lock")
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) {
    return(data.frame(n = n, rep = rep, worker = worker, status = "locked",
                      seconds = 0, path = path, detail = "artifact lock exists",
                      stringsAsFactors = FALSE))
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  warnings_seen <- character(); started <- proc.time()[[3L]]
  answer <- tryCatch(withCallingHandlers(
    run_one_proxlq_replication(n, rep, cfg, parent_root),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  elapsed <- proc.time()[[3L]] - started
  if (inherits(answer, "error")) {
    failure <- list(
      metadata = list(config_id = config_id, n = n, rep = rep, worker = worker,
                      failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
      error = conditionMessage(answer),
      call = paste(deparse(conditionCall(answer)), collapse = " "),
      warnings = unique(warnings_seen), session = capture.output(sessionInfo())
    )
    atomic_save_rds(failure, proxlq_failure_path(path))
    return(data.frame(n = n, rep = rep, worker = worker, status = "failed",
                      seconds = elapsed, path = path, detail = conditionMessage(answer),
                      stringsAsFactors = FALSE))
  }
  answer$warnings <- unique(warnings_seen)
  atomic_save_rds(answer, path)
  if (!validate_proxlq_part(path, cfg, n, rep, parent_root)) {
    return(data.frame(n = n, rep = rep, worker = worker,
                      status = "invalid_after_save", seconds = elapsed, path = path,
                      detail = "saved artifact failed validation", stringsAsFactors = FALSE))
  }
  failure <- proxlq_failure_path(path)
  if (file.exists(failure)) unlink(failure)
  data.frame(n = n, rep = rep, worker = worker, status = "completed",
             seconds = elapsed, path = path, detail = "", stringsAsFactors = FALSE)
}
