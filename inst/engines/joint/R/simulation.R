validate_requested_design <- function(cfg, n = NULL) {
  assert_true(identical(cfg$methods,
                        c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2", "ProxLL-US")),
              "Method list mismatch")
  assert_true(identical(cfg$pdr_methods,
                        c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2")),
              "PDR method list mismatch")
  assert_true(identical(cfg$estimands, c("dose_response", "derivative")),
              "Estimand list mismatch")
  assert_true(identical(cfg$scenarios$scenario, c("CC", "CM", "MC", "MM")),
              "Scenario list mismatch")
  assert_true(cfg$kernel == "epanechnikov", "Only the Epanechnikov kernel is implemented")
  assert_true(cfg$folds >= 2L, "Cross-fitting needs at least two folds")
  assert_true(length(cfg$point_grid) == 9L, "Nine CI reporting points are required")
  assert_true(all(diff(cfg$curve_grid) > 0) && all(diff(cfg$point_grid) > 0),
              "Evaluation grids must increase")
  assert_true(all(grid_membership(cfg$point_grid, cfg$curve_grid)),
              "Point grid must be a subset of the band grid")
  assert_true(max(abs(cfg$point_grid + rev(cfg$point_grid))) < 1e-14,
              "Point grid must be symmetric")
  assert_true(max(abs(range(cfg$curve_grid) - c(-0.9, 0.9))) < 1e-14,
              "Study grid must be [-0.9,0.9]")
  assert_true(identical(cfg$gaussian_band_draws, if (cfg$profile == "production") 10000L
                               else cfg$gaussian_band_draws),
              "Production bands require B=10000")
  if (!is.null(n)) assert_true(n >= 5L * cfg$cv_min_neighbors,
                               "n is too small for configured local fits")
  invisible(TRUE)
}

selection_value <- function(selection, field, scenario) selection[[field]][[scenario]]

bandwidth_table <- function(n, rep, data, cfg, level, derivative,
                            level_plugin, derivative_plugin, h_by_estimand) {
  scale <- robust_exposure_scale(data$A)
  rows <- list(); index <- 0L
  for (estimand in cfg$estimands) for (scenario in cfg$scenarios$scenario) {
    for (method in cfg$methods) {
      index <- index + 1L
      is_level <- estimand == "dose_response"
      selected_h <- h_by_estimand[[estimand]][[method]][[scenario]]
      source_h <- selected_risk <- second_risk <- risk_gap <- normalized_gap <- NA_real_
      lower_hit <- upper_hit <- boundary_hit <- NA
      h_raw <- pilot_h <- integrated_bias2 <- integrated_variance <- NA_real_
      if (method == "PDR-PI1") {
        plugin <- if (is_level) level_plugin[[scenario]] else derivative_plugin[[scenario]]
        h_raw <- plugin$h_raw; pilot_h <- plugin$pilot_h
        integrated_bias2 <- plugin$integrated_bias2
        integrated_variance <- plugin$integrated_variance
        boundary_hit <- plugin$clipped
      } else if (method == "PDR-PI2") {
        source_h <- if (is_level) level$h0 else derivative$h0
      } else if (method == "PDR-CV1") {
        selection <- if (is_level) level$cv1_pdr else derivative$cv1_pdr
        source_h <- selected_h
        selected_risk <- selection_value(selection, "selected_risk", scenario)
        second_risk <- selection_value(selection, "second_risk", scenario)
        risk_gap <- selection_value(selection, "risk_gap", scenario)
        normalized_gap <- selection_value(selection, "normalized_risk_gap", scenario)
        lower_hit <- selection_value(selection, "lower_hit", scenario)
        upper_hit <- selection_value(selection, "upper_hit", scenario)
        boundary_hit <- selection_value(selection, "boundary_hit", scenario)
      } else if (method == "PDR-CV2") {
        selection <- if (is_level) level$cv2_pdr else derivative$cv2_pdr
        source_h <- selected_h
        selected_risk <- selection_value(selection, "selected_risk", scenario)
        second_risk <- selection_value(selection, "second_risk", scenario)
        risk_gap <- selection_value(selection, "risk_gap", scenario)
        normalized_gap <- selection_value(selection, "normalized_risk_gap", scenario)
        lower_hit <- selection_value(selection, "lower_hit", scenario)
        upper_hit <- selection_value(selection, "upper_hit", scenario)
        boundary_hit <- selection_value(selection, "boundary_hit", scenario)
      } else {
        selection <- if (is_level) level$cv1_ll else derivative$cv1_ll
        source_h <- selection$h[[scenario]]
        selected_risk <- selection_value(selection, "selected_risk", scenario)
        second_risk <- selection_value(selection, "second_risk", scenario)
        risk_gap <- selection_value(selection, "risk_gap", scenario)
        normalized_gap <- selection_value(selection, "normalized_risk_gap", scenario)
        lower_hit <- selection_value(selection, "lower_hit", scenario)
        upper_hit <- selection_value(selection, "upper_hit", scenario)
        boundary_hit <- selection_value(selection, "boundary_hit", scenario)
      }
      rate <- if (is_level) -1 / 5 else cfg$derivative_rate_exponent
      rows[[index]] <- data.frame(
        n = n, rep = rep, estimand = estimand, scenario = scenario, method = method,
        h = selected_h, b = if (startsWith(method, "PDR-")) selected_h else NA_real_,
        source_h = source_h, h_raw = h_raw, pilot_h = pilot_h,
        integrated_bias2 = integrated_bias2, integrated_variance = integrated_variance,
        selected_risk = selected_risk, second_risk = second_risk,
        risk_gap = risk_gap, normalized_risk_gap = normalized_gap,
        lower_hit = lower_hit, upper_hit = upper_hit, boundary_hit = boundary_hit,
        rate_exponent = rate, normalized_h = selected_h / (scale * n^rate),
        n_h3 = n * selected_h^3, n_h7 = n * selected_h^7,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

candidate_allowed <- function(h, allowed) {
  vapply(h, function(value) any(abs(log(value / allowed)) < 1e-12), logical(1L))
}

cv_diagnostic_table <- function(n, rep, cfg, level, derivative) {
  rows <- list(); index <- 0L
  for (estimand in cfg$estimands) {
    object <- if (estimand == "dose_response") level else derivative
    bandwidth <- object$final$bandwidth
    for (scenario in cfg$scenarios$scenario) {
      if (estimand == "dose_response") {
        risk_pdr <- object$final$risk_pdr[, scenario]
        risk_prox <- object$final$risk_ll[, scenario]
        invalid_pdr <- invalid_prox <- object$final$invalid_evaluations
        evaluations <- object$final$trimmed_evaluations
        min_neighbors <- min_pivot <- rep(NA_real_, length(bandwidth))
        density_floor_fraction <- 0
      } else {
        risk_pdr <- object$final$risk_pdr[, scenario]
        risk_prox <- object$final$risk_ll[, scenario]
        invalid_pdr <- object$final$invalid_pdr
        invalid_prox <- object$final$invalid_ll
        evaluations <- object$final$tapered_evaluations
        min_neighbors <- object$final$min_neighbors
        min_pivot <- object$final$min_pivot_ratio
        density_floor_fraction <- object$density$floor_fraction
      }
      for (selector in c("CV1-PDR", "CV2-PDR", "CV1-Prox")) {
        index <- index + 1L
        if (selector == "CV1-PDR") {
          selection <- object$cv1_pdr
          risk <- risk_pdr; invalid <- invalid_pdr
        } else if (selector == "CV2-PDR") {
          selection <- object$cv2_pdr
          risk <- risk_pdr; invalid <- invalid_pdr
        } else {
          selection <- object$cv1_ll; risk <- risk_prox; invalid <- invalid_prox
        }
        candidate_rows <- data.frame(
          n = n, rep = rep, estimand = estimand, scenario = scenario,
          selector = selector, bandwidth = bandwidth, risk = risk,
          allowed = candidate_allowed(bandwidth, selection$allowed[[scenario]]),
          finite_candidate = is.finite(risk),
          selected = abs(log(bandwidth / selection$h[[scenario]])) < 1e-12,
          invalid_evaluations = invalid, target_evaluations = evaluations,
          min_neighbors = min_neighbors, min_pivot_ratio = min_pivot,
          density_floor_fraction = density_floor_fraction,
          stringsAsFactors = FALSE
        )
        rows[[index]] <- candidate_rows[candidate_rows$allowed, , drop = FALSE]
      }
    }
  }
  do.call(rbind, rows)
}

run_one_replication <- function(n, rep, cfg, root = project_root()) {
  validate_requested_design(cfg, n)
  started <- proc.time()[[3L]]
  dgp_seed <- seed_for(n, rep, 1L, cfg$base_seed)
  fold_seed <- seed_for(n, rep, 2L, cfg$base_seed)
  data <- simulate_proximal_data(n, dgp_seed)
  fold <- balanced_folds(n, cfg$folds, fold_seed)

  nuisance_started <- proc.time()[[3L]]
  crossfit <- crossfit_bridges(data, fold, cfg)
  nuisance_seconds <- proc.time()[[3L]] - nuisance_started

  bandwidth_started <- proc.time()[[3L]]
  level <- select_level_bandwidths(data$A, crossfit$xi, cfg)
  derivative <- select_derivative_bandwidths(data$A, crossfit$xi, cfg)
  level_bounds <- range(level$initial); derivative_bounds <- range(derivative$initial)
  level_plugin <- lapply(colnames(crossfit$xi), function(scenario)
    level_plugin_bandwidth(data$A, crossfit$xi[, scenario],
                           level$cv1_ll$h[[scenario]], level_bounds, cfg))
  derivative_plugin <- lapply(colnames(crossfit$xi), function(scenario)
    derivative_plugin_bandwidth(data$A, crossfit$xi[, scenario], derivative_bounds, cfg))
  names(level_plugin) <- names(derivative_plugin) <- colnames(crossfit$xi)
  level_pi1 <- vapply(level_plugin, `[[`, numeric(1L), "h")
  derivative_pi1 <- vapply(derivative_plugin, `[[`, numeric(1L), "h")
  scenarios <- cfg$scenarios$scenario
  h_by_estimand <- list(
    dose_response = list(
      `PDR-PI1` = level_pi1,
      `PDR-PI2` = setNames(rep(level$h0, length(scenarios)), scenarios),
      `PDR-CV1` = level$cv1_pdr$h,
      `PDR-CV2` = level$cv2_pdr$h,
      `ProxLL-US` = level$cv1_ll$h / log10(n)
    ),
    derivative = list(
      `PDR-PI1` = derivative_pi1,
      `PDR-PI2` = setNames(rep(derivative$h0, length(scenarios)), scenarios),
      `PDR-CV1` = derivative$cv1_pdr$h,
      `PDR-CV2` = derivative$cv2_pdr$h,
      `ProxLL-US` = derivative$cv1_ll$h / log10(n)
    )
  )
  bandwidth_seconds <- proc.time()[[3L]] - bandwidth_started

  inference_started <- proc.time()[[3L]]
  curves <- bands <- list(); curve_index <- band_index <- 0L
  for (scenario_index in seq_along(scenarios)) {
    scenario <- scenarios[[scenario_index]]
    gaussian <- shared_gaussian_draws(
      cfg$gaussian_band_draws, length(cfg$curve_grid),
      seed_for(n, rep, 300L + scenario_index, cfg$base_seed)
    )
    for (estimand in cfg$estimands) for (method in cfg$methods) {
      pointwise <- evaluate_pointwise(
        data, crossfit, scenario, method, estimand,
        h_by_estimand[[estimand]][[method]][[scenario]], cfg
      )
      curve_index <- curve_index + 1L
      if (method %in% cfg$pdr_methods) {
        answer <- attach_gaussian_band(pointwise, cfg, gaussian)
        curves[[curve_index]] <- answer$curve
        band_index <- band_index + 1L; bands[[band_index]] <- answer$band
      } else curves[[curve_index]] <- pointwise$curve
    }
    rm(gaussian)
  }
  inference_seconds <- proc.time()[[3L]] - inference_started
  curve <- do.call(rbind, curves)
  curve$n <- n; curve$rep <- rep
  curve$pointwise_coverage <- curve$truth >= curve$ci_lower & curve$truth <= curve$ci_upper
  curve$ci_width <- curve$ci_upper - curve$ci_lower
  curve$squared_error <- (curve$estimate - curve$truth)^2
  curve <- curve[c("n", "rep", setdiff(names(curve), c("n", "rep")))]
  rownames(curve) <- NULL
  band <- do.call(rbind, bands)
  band$n <- n; band$rep <- rep
  band <- band[c("n", "rep", setdiff(names(band), c("n", "rep")))]
  rownames(band) <- NULL

  bandwidth <- bandwidth_table(n, rep, data, cfg, level, derivative,
                               level_plugin, derivative_plugin, h_by_estimand)
  cv_diagnostics <- cv_diagnostic_table(n, rep, cfg, level, derivative)
  list(
    metadata = list(
      version = cfg$version, profile = cfg$profile, config_id = config_fingerprint(cfg),
      n = n, rep = rep, dgp_seed = dgp_seed, fold_seed = fold_seed,
      gaussian_band_draws = cfg$gaussian_band_draws,
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      session = capture.output(sessionInfo())
    ),
    curve = curve, band = band, bandwidth = bandwidth,
    nuisance_diagnostics = transform(crossfit$diagnostics, n = n, rep = rep),
    cv_diagnostics = cv_diagnostics,
    dgp_diagnostics = dgp_diagnostics(data, range(cfg$curve_grid)),
    timing = c(nuisance_seconds = nuisance_seconds,
               bandwidth_seconds = bandwidth_seconds,
               inference_seconds = inference_seconds,
               total_seconds = proc.time()[[3L]] - started)
  )
}

result_path <- function(outdir, config_id, n, rep) {
  file.path(outdir, paste0("config_", config_id), paste0("n_", n),
            sprintf("rep_%04d.rds", rep))
}

validate_part <- function(path, cfg, n, rep) {
  if (!file.exists(path)) return(FALSE)
  answer <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(answer) || is.null(answer$metadata)) return(FALSE)
  metadata_ok <- identical(answer$metadata$config_id, config_fingerprint(cfg)) &&
    identical(as.integer(answer$metadata$n), as.integer(n)) &&
    identical(as.integer(answer$metadata$rep), as.integer(rep)) &&
    identical(as.integer(answer$metadata$gaussian_band_draws),
              as.integer(cfg$gaussian_band_draws))
  if (!metadata_ok || !is.data.frame(answer$curve) || !is.data.frame(answer$band) ||
      !is.data.frame(answer$bandwidth) || !is.data.frame(answer$nuisance_diagnostics) ||
      !is.data.frame(answer$cv_diagnostics)) return(FALSE)
  curve <- answer$curve
  expected_rows <- length(cfg$estimands) * nrow(cfg$scenarios) * length(cfg$methods) *
    length(cfg$curve_grid)
  key <- paste(curve$estimand, curve$scenario, curve$method,
               sprintf("%.12g", curve$a), sep = "|")
  curve_ok <- nrow(curve) == expected_rows && !anyDuplicated(key) &&
    all(curve$n == n) && all(curve$rep == rep) &&
    all(is.finite(as.matrix(curve[c("truth", "estimate", "se", "ci_lower", "ci_upper",
                                      "bandwidth", "phi_mean", "max_condition",
                                      "reproduction_error", "local_n", "ci_width",
                                      "squared_error")]))) &&
    all(curve$se > 0) && all(curve$local_n >= cfg$inference_min_neighbors) &&
    max(curve$reproduction_error) <= cfg$local_reproduction_tolerance &&
    all(is.na(curve$band_lower[curve$method == "ProxLL-US"])) &&
    all(is.na(curve$band_upper[curve$method == "ProxLL-US"])) &&
    all(is.finite(curve$band_lower[curve$method != "ProxLL-US"])) &&
    all(is.finite(curve$band_upper[curve$method != "ProxLL-US"]))
  band <- answer$band
  band_key <- paste(band$estimand, band$scenario, band$method, sep = "|")
  band_ok <- nrow(band) == length(cfg$estimands) * nrow(cfg$scenarios) *
    length(cfg$pdr_methods) && !anyDuplicated(band_key) &&
    all(is.finite(as.matrix(band[c("critical", "mean_width", "min_eigen",
                                     "projection_norm")])))
  bandwidth <- answer$bandwidth
  bandwidth_key <- paste(bandwidth$estimand, bandwidth$scenario, bandwidth$method, sep = "|")
  bandwidth_ok <- nrow(bandwidth) == length(cfg$estimands) * nrow(cfg$scenarios) *
    length(cfg$methods) && !anyDuplicated(bandwidth_key) &&
    all(is.finite(bandwidth$h)) && all(bandwidth$h > 0) &&
    all(abs(bandwidth$b[!is.na(bandwidth$b)] - bandwidth$h[!is.na(bandwidth$b)]) < 1e-14) &&
    all(is.finite(as.matrix(bandwidth[c("normalized_h", "n_h3", "n_h7")])))
  nuisance_ok <- nrow(answer$nuisance_diagnostics) == cfg$folds * 4L &&
    all(answer$nuisance_diagnostics$converged)
  cv <- answer$cv_diagnostics
  cv_group <- interaction(cv$estimand, cv$scenario, cv$selector, drop = TRUE)
  expected_cv_groups <- length(cfg$estimands) * nrow(cfg$scenarios) * 3L
  cv_ok <- nrow(cv) >= expected_cv_groups &&
    length(unique(cv_group)) == expected_cv_groups &&
    all(cv$allowed) &&
    all(vapply(split(cv$selected, cv_group), function(value) sum(value) == 1L,
               logical(1L))) &&
    all(cv$finite_candidate[cv$selected]) &&
    all(is.finite(cv$bandwidth)) && all(cv$bandwidth > 0)
  timing_ok <- length(answer$timing) == 4L && all(is.finite(answer$timing)) &&
    all(answer$timing >= 0)
  isTRUE(curve_ok && band_ok && bandwidth_ok && nuisance_ok && cv_ok && timing_ok)
}
