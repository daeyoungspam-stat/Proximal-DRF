simulation_conditions <- function(cfg) {
  grid <- expand.grid(degree = cfg$null_degrees, delta = cfg$delta_grid,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[order(grid$degree, grid$delta), , drop = FALSE]
  grid$null <- unname(cfg$null_labels[as.character(grid$degree)])
  grid$condition <- sprintf("d%d_delta_%s", grid$degree,
                            gsub("[.]", "p", sprintf("%.6f", grid$delta)))
  grid[c("condition", "null", "degree", "delta")]
}

validate_requested_design <- function(cfg, n = NULL) {
  assert_true(identical(cfg$methods,
                        c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                          "ProxLL-raw", "ProxLL-raw-matched",
                          "ProxLL-US", "ProxLL-US-matched")),
              "Inference-method list mismatch")
  assert_true(identical(cfg$bandwidth_methods,
                        c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                          "ProxLL-raw", "ProxLL-US")),
              "Bandwidth-method list mismatch")
  assert_true(identical(cfg$null_degrees, 0:2), "Null-degree mismatch")
  assert_true(identical(cfg$multiplier_laws, "two_point"),
              "Multiplier-law mismatch")
  assert_true(identical(cfg$bootstrap_schemes,
                        c("fold_parametric", "fold_nonparametric")),
              "Bootstrap-scheme mismatch")
  expected_scenarios <- if (cfg$profile == "calibration")
    c("CC", "CM", "MC") else c("CC", "CM", "MC", "MM")
  assert_true(identical(cfg$scenarios$scenario, expected_scenarios),
              "Nuisance-scenario mismatch")
  assert_true(cfg$bootstrap_draws >= 19L, "Too few bootstrap draws")
  assert_true(abs(sum(cfg$quadrature_weight) - 1.8) < 1e-12,
              "Quadrature does not integrate constants on [-0.9,0.9]")
  assert_true(identical(cfg$quadrature_order, 512L) &&
                identical(cfg$quadrature_panels, 64L) &&
                identical(cfg$quadrature_nodes_per_panel, 8L) &&
                length(cfg$quadrature_node) == cfg$quadrature_order,
              "Production integration rule is not the frozen composite rule")
  assert_true(all(cfg$weighted_quadrature >= 0), "Negative weighted quadrature")
  if (!is.null(n)) assert_true(n >= 5L * cfg$cv_min_neighbors,
                               "n is too small for configured local fits")
  invisible(TRUE)
}

bandwidth_map <- function(column, level, plugin, n, cfg) {
  answer <- c(
    `PDR-PI1` = plugin[[column]]$h,
    `PDR-PI2` = level$h0,
    `PDR-CV1` = level$cv1_pdr$h[[column]],
    `PDR-CV2` = level$cv2_pdr$h[[column]],
    `ProxLL-raw` = level$cv1_ll$h[[column]],
    `ProxLL-US` = level$cv1_ll$h[[column]] / log10(n)
  )
  answer[cfg$bandwidth_methods]
}

selection_diagnostics <- function(column, level, plugin, h, key, condition_row) {
  cv1_pdr_allowed <- level$cv1_pdr$allowed[[column]]
  cv2_pdr_allowed <- level$cv2_pdr$allowed[[column]]
  cv1_ll_allowed <- level$cv1_ll$allowed[[column]]
  pi_bounds <- range(level$initial)
  data.frame(
    condition = key$condition, null = condition_row$null,
    degree = condition_row$degree, delta = condition_row$delta,
    scenario = key$scenario,
    method = names(h), bandwidth = unname(h),
    source_bandwidth = c(plugin[[column]]$h_raw, level$h0,
                         level$cv1_pdr$h[[column]], level$cv2_pdr$h[[column]],
                         level$cv1_ll$h[[column]], level$cv1_ll$h[[column]]),
    used_to_source_ratio = unname(h) / c(plugin[[column]]$h_raw, level$h0,
                         level$cv1_pdr$h[[column]], level$cv2_pdr$h[[column]],
                         level$cv1_ll$h[[column]], level$cv1_ll$h[[column]]),
    selector = c("PI1", "PI2-fixed", "CV1-PDR", "CV2-PDR",
                 "CV1-ProxLL", "CV1-ProxLL-source"),
    boundary_applicable = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    lower_hit = c(plugin[[column]]$lower_hit, NA, level$cv1_pdr$lower_hit[[column]],
                  level$cv2_pdr$lower_hit[[column]],
                  level$cv1_ll$lower_hit[[column]], level$cv1_ll$lower_hit[[column]]),
    upper_hit = c(plugin[[column]]$upper_hit, NA, level$cv1_pdr$upper_hit[[column]],
                  level$cv2_pdr$upper_hit[[column]],
                  level$cv1_ll$upper_hit[[column]], level$cv1_ll$upper_hit[[column]]),
    selected_risk = c(NA, NA, level$cv1_pdr$selected_risk[[column]],
                      level$cv2_pdr$selected_risk[[column]],
                      level$cv1_ll$selected_risk[[column]],
                      level$cv1_ll$selected_risk[[column]]),
    second_risk = c(NA, NA, level$cv1_pdr$second_risk[[column]],
                    level$cv2_pdr$second_risk[[column]],
                    level$cv1_ll$second_risk[[column]],
                    level$cv1_ll$second_risk[[column]]),
    risk_gap = c(NA, NA, level$cv1_pdr$risk_gap[[column]],
                 level$cv2_pdr$risk_gap[[column]],
                 level$cv1_ll$risk_gap[[column]], level$cv1_ll$risk_gap[[column]]),
    normalized_risk_gap = c(NA, NA, level$cv1_pdr$normalized_risk_gap[[column]],
                            level$cv2_pdr$normalized_risk_gap[[column]],
                            level$cv1_ll$normalized_risk_gap[[column]],
                            level$cv1_ll$normalized_risk_gap[[column]]),
    selector_min = c(pi_bounds[[1L]], NA, min(cv1_pdr_allowed), min(cv2_pdr_allowed),
                     min(cv1_ll_allowed), min(cv1_ll_allowed)),
    selector_max = c(pi_bounds[[2L]], NA, max(cv1_pdr_allowed), max(cv2_pdr_allowed),
                     max(cv1_ll_allowed), max(cv1_ll_allowed)),
    selector_count = c(NA_integer_, 1L, length(cv1_pdr_allowed), length(cv2_pdr_allowed),
                       length(cv1_ll_allowed), length(cv1_ll_allowed)),
    stringsAsFactors = FALSE
  )
}

run_one_replication <- function(n, rep, cfg, root = project_root()) {
  context <- checkpoint_context(n, rep, cfg, root)
  checkpoint_check_deadline(context)
  withCallingHandlers(
    run_one_replication_impl(n, rep, cfg, root, context),
    warning = function(w) {
      context$warnings <- unique(c(context$warnings, conditionMessage(w)))
    })
}

run_one_replication_impl <- function(n, rep, cfg, root, context) {
  validate_requested_design(cfg, n); started <- proc.time()[[3L]]
  conditions <- simulation_conditions(cfg)
  dgp_seed <- seed_for(n, rep, 1L, cfg$base_seed)
  fold_seed <- seed_for(n, rep, 2L, cfg$base_seed)
  multiplier_seed <- seed_for(n, rep, 500L, cfg$base_seed)
  base <- simulate_proximal_base(n, dgp_seed, keep_latent = TRUE)
  outcome <- outcome_matrix(base, conditions, cfg)
  fold <- balanced_folds(n, cfg$folds, fold_seed)

  nuisance_started <- proc.time()[[3L]]
  crossfit <- crossfit_all_conditions(base, outcome, fold, cfg)
  nuisance_seconds <- proc.time()[[3L]] - nuisance_started

  bandwidth_started <- proc.time()[[3L]]
  level <- select_level_bandwidths(base$A, crossfit$xi, cfg)
  original_level <- level
  feasibility <- NULL
  if (isTRUE(cfg$fold_cv_filter)) {
    feasibility <- new_fold_feasibility(base$A, fold, cfg)
    level <- filter_level_cv(level, base$A, fold, cfg, feasibility)
  }
  bounds <- range(level$initial)
  plugin <- lapply(seq_len(ncol(crossfit$xi)), function(column) {
    level_plugin_bandwidth(base$A, crossfit$xi[, column],
                           original_level$cv1_ll$h[[column]], bounds, cfg)
  })
  h_by_column <- lapply(seq_len(ncol(crossfit$xi)), bandwidth_map,
                        level = level, plugin = plugin, n = n, cfg = cfg)
  bandwidth_seconds <- proc.time()[[3L]] - bandwidth_started

  if (isTRUE(cfg$fold_cv_filter)) {
    original_h <- lapply(seq_len(ncol(crossfit$xi)), bandwidth_map,
      level = original_level, plugin = plugin, n = n, cfg = cfg)
    failed <- list()
    for (column in seq_along(h_by_column)) for (method in cfg$bandwidth_methods) {
      degree <- if (method %in% cfg$pdr_methods) 2L else 1L
      check <- feasibility$check(h_by_column[[column]][[method]], degree)
      if (!check$feasible) failed[[length(failed) + 1L]] <- list(
        column = column, key = crossfit$key[column, ], method = method, check = check)
    }
    audit <- list(n = n, rep = rep, key = crossfit$key,
      original_h = original_h, filtered_h = h_by_column,
      selection = level$admissibility, checks = as.list(feasibility$cache),
      remaining_failures = failed,
      dgp_seed = dgp_seed, fold_seed = fold_seed, multiplier_seed = multiplier_seed,
      preparation_seconds = proc.time()[[3L]] - started)
    audit_dir <- Sys.getenv("PDR_BANDWIDTH_AUDIT_DIR")
    if (nzchar(audit_dir)) atomic_save_rds(audit,
      file.path(audit_dir, sprintf("rep_%04d_bandwidth_audit.rds", rep)), compress = FALSE)
    if (length(failed)) stop("Non-CV or remaining operator failure after CV filtering; see bandwidth audit")
    if (identical(Sys.getenv("PDR_PREP_ONLY"), "1")) return(audit)
  }

  checkpoint_prepare(context, crossfit$key, list(
    base = base, outcome = outcome, fold = fold, xi = crossfit$xi,
    nuisance = crossfit$diagnostics, bandwidth = h_by_column), cfg)

  inference_started <- proc.time()[[3L]]
  random_blocks <- make_random_index_blocks(n, cfg, multiplier_seed)
  validate_random_blocks(random_blocks, n, cfg)
  fold_indices <- lapply(seq_len(cfg$folds), function(k) which(fold == k))
  fold_projection_by_degree <- lapply(fold_indices, function(index)
    setNames(lapply(cfg$null_degrees, function(degree)
      null_projection(base$A[index], degree)), as.character(cfg$null_degrees)))
  request_keys <- unlist(lapply(h_by_column, function(h) {
    vapply(cfg$bandwidth_methods, function(method) {
      degree <- if (method %in% cfg$pdr_methods) 2L else 1L
      operator_cache_key(degree, h[[method]])
    }, character(1L))
  }), use.names = FALSE)
  key_count <- table(request_keys)
  cacheable_keys <- names(key_count)[key_count > 1L]
  operator_cache <- lapply(seq_len(cfg$folds), function(k)
    new.env(parent = emptyenv()))
  results <- bandwidth_rows <- residual_rows <- vector("list", ncol(crossfit$xi))
  if (!is.null(context$saved)) {
    results <- context$saved$results
    bandwidth_rows <- context$saved$bandwidth_rows
    residual_rows <- context$saved$residual_rows
  }
  current_timing <- function() c(nuisance = nuisance_seconds,
    bandwidth = bandwidth_seconds, inference = proc.time()[[3L]] - inference_started,
    total = proc.time()[[3L]] - started)
  # Commit the regenerated preparation before a deadline pause, including the
  # compute time and any new warnings from this segment (queue time excluded).
  checkpoint_commit(context, context$completed, results, bandwidth_rows,
                    residual_rows, current_timing())
  checkpoint_check_deadline(context)
  for (column in seq_len(ncol(crossfit$xi))) {
    if (column <= context$completed) next
    key <- crossfit$key[column, , drop = FALSE]
    condition_row <- conditions[match(key$condition, conditions$condition), , drop = FALSE]
    h <- h_by_column[[column]]
    fold_operators <- lapply(seq_len(cfg$folds), function(k) {
      index <- fold_indices[[k]]
      build_operator_set(base$A[index], h, cfg, cache = operator_cache[[k]],
                         cacheable_keys = cacheable_keys,
                         min_neighbors = cfg$fold_inference_min_neighbors)
    })
    # At observed doses, residualize against the estimator for this method;
    # unlike the null fit, this fit varies with degree and bandwidth.
    np_fits <- lapply(seq_len(cfg$folds), function(k) {
      index <- fold_indices[[k]]
      values <- lapply(cfg$bandwidth_methods, function(method) {
        polynomial_degree <- if (method %in% cfg$pdr_methods) 2L else 1L
        nonparametric_fitted_values(
          base$A[index], crossfit$xi[index, column], h[[method]],
          polynomial_degree, cfg)
      })
      setNames(values, cfg$bandwidth_methods)
    })
    degree_key <- as.character(condition_row$degree[[1L]])
    test <- wild_bootstrap_tests(
      y = crossfit$xi[, column], a = base$A,
      degree = condition_row$degree[[1L]], cfg = cfg, fold = fold,
      fold_operators = fold_operators, np_fits = np_fits,
      random_blocks = random_blocks,
      random_blocks_validated = TRUE,
      fold_projections = lapply(fold_projection_by_degree, `[[`, degree_key)
    )
    result <- test$result
    result$condition <- key$condition; result$null <- condition_row$null[[1L]]
    result$degree <- condition_row$degree[[1L]]; result$delta <- condition_row$delta[[1L]]
    result$scenario <- key$scenario; result$identified <-
      cfg$scenarios$identified[match(key$scenario, cfg$scenarios$scenario)]
    # The engine records the actual scheme- and method-specific residual
    # diagnostics on each test row; do not overwrite them with null residuals.
    results[[column]] <- result
    bandwidth_rows[[column]] <- selection_diagnostics(column, level, plugin, h, key,
                                                       condition_row)
    residual_rows[[column]] <- data.frame(
      condition = key$condition, null = condition_row$null,
      degree = condition_row$degree, delta = condition_row$delta,
      scenario = key$scenario,
      as.list(test$residual_diagnostics), stringsAsFactors = FALSE
    )
    checkpoint_commit(context, column, results, bandwidth_rows, residual_rows,
                      current_timing())
    if (identical(Sys.getenv("PDR_PROGRESS"), "1") &&
        (column %% 12L == 0L || column == ncol(crossfit$xi))) {
      message(sprintf("n=%d rep=%d conditions=%d/%d elapsed=%.1fs",
                      n, rep, column, ncol(crossfit$xi),
                      proc.time()[[3L]] - started))
    }
    # A condition is indivisible; a hard timeout during one condition loses
    # only that condition plus the uncheckpointed part of preparation.
    checkpoint_check_deadline(context)
  }
  inference_seconds <- proc.time()[[3L]] - inference_started
  accumulated <- context$prior_timing + current_timing()
  test_result <- do.call(rbind, results); test_result$n <- n; test_result$rep <- rep
  test_result <- test_result[c("n", "rep", setdiff(names(test_result), c("n", "rep")))]
  bandwidth <- do.call(rbind, bandwidth_rows); bandwidth$n <- n; bandwidth$rep <- rep
  bandwidth <- bandwidth[c("n", "rep", setdiff(names(bandwidth), c("n", "rep")))]
  residual <- do.call(rbind, residual_rows); residual$n <- n; residual$rep <- rep
  residual <- residual[c("n", "rep", setdiff(names(residual), c("n", "rep")))]
  list(
    metadata = list(version = cfg$version, profile = cfg$profile,
                    config_id = context$config_id, n = n, rep = rep,
                    dgp_seed = dgp_seed, fold_seed = fold_seed,
                    multiplier_seed = multiplier_seed,
                    bootstrap_draws = cfg$bootstrap_draws,
                    multiplier_laws = cfg$multiplier_laws,
                    bootstrap_schemes = cfg$bootstrap_schemes,
                    folds = cfg$folds, fold_sizes = lengths(fold_indices),
                    bandwidth_selection = if (isTRUE(cfg$fold_cv_filter))
                      cfg$fold_cv_filter_rule else "pooled_original_level_selectors_fixed_in_all_folds_and_draws",
                    nonparametric_residual_fit = cfg$nonparametric_residual_fit,
                    singular_observation_fit = cfg$singular_observation_fit,
                    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                    checkpoint_segments = context$segments,
                    session = capture.output(sessionInfo())),
    warnings = unique(context$warnings),
    test = test_result, bandwidth = bandwidth, residual_diagnostics = residual,
    nuisance_diagnostics = transform(crossfit$diagnostics, n = n, rep = rep),
    dgp_diagnostics = data.frame(n = n, rep = rep,
                                 as.list(dgp_diagnostics(base, cfg)),
                                 check.names = FALSE),
    timing = data.frame(n = n, rep = rep,
                        nuisance_seconds = unname(accumulated[["nuisance"]]),
                        bandwidth_seconds = unname(accumulated[["bandwidth"]]),
                        inference_seconds = unname(accumulated[["inference"]]),
                        total_seconds = unname(accumulated[["total"]]))
  )
}
