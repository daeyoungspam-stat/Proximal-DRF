make_random_index_blocks <- function(n, cfg, multiplier_seed) {
  with_preserved_seed(multiplier_seed, {
    answer <- list(); start <- 1L
    upper <- 2147483647
    while (start <= cfg$bootstrap_draws) {
      stop_at <- min(cfg$bootstrap_draws, start + cfg$bootstrap_block - 1L)
      count <- stop_at - start + 1L
      value <- floor(runif(n * count, min = 0, max = upper + 1))
      answer[[length(answer) + 1L]] <- matrix(as.integer(value), nrow = n)
      start <- stop_at + 1L
    }
    answer
  })
}

multiplier_from_index <- function(index, law = "two_point") {
  assert_true(identical(law, "two_point"), "Only two-point multipliers are configured")
  u <- (as.numeric(index) + 0.5) / 2147483648
  dim(u) <- dim(index)
  root5 <- sqrt(5)
  ifelse(u < (5 + root5) / 10, (1 - root5) / 2, (1 + root5) / 2)
}

make_multiplier_blocks <- function(n, cfg, multiplier_seed, law = "two_point") {
  lapply(make_random_index_blocks(n, cfg, multiplier_seed),
         multiplier_from_index, law = law)
}

validate_random_blocks <- function(blocks, n, cfg) {
  assert_true(length(blocks) > 0L &&
                sum(vapply(blocks, ncol, integer(1L))) == cfg$bootstrap_draws,
              "Cached random blocks have the wrong number of draws")
  assert_true(all(vapply(blocks, nrow, integer(1L)) == n),
              "Cached random blocks have the wrong sample size")
  assert_true(all(vapply(blocks, function(x)
    is.integer(x) && !anyNA(x) && all(x >= 0L), logical(1L))),
    "Cached random blocks have an invalid integer representation")
  invisible(TRUE)
}

residual_diagnostics <- function(residual) {
  centered <- residual - mean(residual)
  variance <- mean(centered^2)
  excess_kurtosis <- if (variance > 0)
    mean(centered^4) / variance^2 - 3 else NA_real_
  c(residual_sd = sd(residual),
    residual_abs_q99 = unname(quantile(abs(residual), .99, type = 8)),
    residual_abs_max = max(abs(residual)),
    residual_excess_kurtosis = excess_kurtosis)
}

fold_bootstrap_finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

fold_average_bootstrap_state <- function(y, a, degree, cfg, fold,
                                         fold_operators,
                                         fold_projections = NULL,
                                         np_fits = NULL) {
  n <- length(y)
  assert_true(length(a) == n && length(fold) == n &&
                all(is.finite(y)) && all(is.finite(a)),
              "Invalid observations supplied to fold bootstrap")
  assert_true(identical(sort(unique(as.integer(fold))), seq_len(cfg$folds)),
              "Cross-fitting fold labels mismatch")
  assert_true(length(fold_operators) == cfg$folds,
              "Fold operator count mismatch")
  if (is.null(fold_projections)) fold_projections <- vector("list", cfg$folds)
  assert_true(length(fold_projections) == cfg$folds,
              "Fold projection count mismatch")
  need_np <- "fold_nonparametric" %in% cfg$bootstrap_schemes
  apply_backend <- cfg$bootstrap_apply_backend %||% "matrix"
  assert_true(apply_backend %in% c("matrix", "blas"),
              "Unknown bootstrap smoother-product backend")
  if (apply_backend != "matrix")
    assert_true(exists("fast_fold_operator_apply", mode = "function", inherits = TRUE) &&
                  exists("blas_dense_apply_cpp", mode = "function", inherits = TRUE),
                "Fast smoother backend must be loaded before bootstrap")
  if (need_np) assert_true(length(np_fits) == cfg$folds,
                           "Nonparametric fold fits are required")
  specs <- inference_specs(degree, cfg)
  pieces <- lapply(seq_len(cfg$folds), function(k) {
    index <- which(fold == k)
    operators <- fold_operators[[k]]
    assert_true(identical(names(operators), cfg$bandwidth_methods),
                "Fold operator order mismatch")
    observed <- observed_test_components(
      y[index], a[index], degree, operators, cfg,
      projection = fold_projections[[k]], specs = specs)
    np_residual <- np_diagnostics <- NULL
    if (need_np) {
      assert_true(identical(names(np_fits[[k]]), cfg$bandwidth_methods),
                  "Nonparametric-fit method order mismatch")
      np_residual <- lapply(np_fits[[k]], function(fit) {
        fitted <- if (is.list(fit)) fit$fitted else fit
        assert_true(length(fitted) == length(index) && all(is.finite(fitted)),
                    "Invalid nonparametric residual fitted values")
        as.numeric(y[index]) - as.numeric(fitted)
      })
      np_diagnostics <- lapply(names(np_residual), function(method) {
        fit <- np_fits[[k]][[method]]
        zero_index <- if (is.list(fit)) fit$zero_index %||% integer() else integer()
        zero_count <- if (is.list(fit)) fit$zero_count %||% length(zero_index) else 0L
        assert_true(length(zero_index) == zero_count &&
                      !anyDuplicated(zero_index) && !anyNA(zero_index) &&
                      all(zero_index >= 1L & zero_index <= length(index)),
                    "Nonparametric zero-convention diagnostics are inconsistent")
        zero_a <- a[index][zero_index]
        list(zero_count = zero_count,
             zero_inside_test_support = sum(zero_a >= cfg$test_support[[1L]] &
                                              zero_a <= cfg$test_support[[2L]]),
             max_original_condition = if (is.list(fit)) fit$max_condition %||% NA_real_ else NA_real_,
             max_solve_condition = if (is.list(fit)) fit$max_solve_condition %||% NA_real_ else NA_real_,
             exact_zero_residual_count = sum(np_residual[[method]] == 0),
             residual_abs_max = max(abs(np_residual[[method]])))
      })
      names(np_diagnostics) <- names(np_residual)
    }
    dense_operators <- if (apply_backend == "blas")
      lapply(operators, function(op) as.matrix(op$S)) else NULL
    list(index = index, n = length(index), weight = length(index) / n,
         observed = observed, operators = operators,
         stacked = if (apply_backend == "matrix")
           as(do.call(rbind, lapply(operators, `[[`, "S")), "dgCMatrix") else NULL,
         dense_operators = dense_operators, apply_backend = apply_backend,
         np_residual = np_residual, np_diagnostics = np_diagnostics)
  })
  weights <- vapply(pieces, `[[`, numeric(1L), "weight")
  rows <- lapply(specs$method, function(method) {
    components <- lapply(pieces, function(p) p$observed$method[[method]])
    op_method <- components[[1L]]$operator_method
    operators <- lapply(pieces, function(p) p$operators[[op_method]])
    h <- vapply(operators, `[[`, numeric(1L), "h")
    assert_true(all(h == h[[1L]]), "Bandwidth must be shared across folds")
    data.frame(
      method = method, operator_method = op_method,
      null_fit = components[[1L]]$null_fit,
      statistic = sum(weights * vapply(components, `[[`, numeric(1L), "statistic")),
      bandwidth = h[[1L]], polynomial_degree = operators[[1L]]$degree,
      null_fit_identity_error = max(vapply(components, `[[`, numeric(1L), "identity_error")),
      null_smoothing_bias_max = max(vapply(components, function(z)
        max(abs(z$smoothing_bias)), numeric(1L))),
      null_smoothing_bias_l2 = sqrt(sum(weights * vapply(components, function(z)
        weighted_integral(z$smoothing_bias^2, cfg), numeric(1L)))),
      max_condition = max(unlist(lapply(operators, `[[`, "condition"))),
      min_local_n = min(unlist(lapply(operators, `[[`, "local_n"))),
      max_reproduction_error = max(unlist(lapply(operators, `[[`, "reproduction_error"))),
      np_zero_count = if (need_np) sum(vapply(pieces, function(p)
        as.numeric(p$np_diagnostics[[op_method]]$zero_count), numeric(1L))) else NA_real_,
      np_zero_inside_test_support = if (need_np) sum(vapply(pieces, function(p)
        as.numeric(p$np_diagnostics[[op_method]]$zero_inside_test_support), numeric(1L))) else NA_real_,
      np_max_original_condition = if (need_np) fold_bootstrap_finite_max(vapply(pieces, function(p)
        p$np_diagnostics[[op_method]]$max_original_condition, numeric(1L))) else NA_real_,
      np_infinite_original_condition_folds = if (need_np) sum(vapply(pieces, function(p)
        is.infinite(p$np_diagnostics[[op_method]]$max_original_condition), logical(1L))) else NA_integer_,
      np_max_solve_condition = if (need_np) fold_bootstrap_finite_max(vapply(pieces, function(p)
        p$np_diagnostics[[op_method]]$max_solve_condition, numeric(1L))) else NA_real_,
      np_exact_zero_residual_count = if (need_np) sum(vapply(pieces, function(p)
        as.numeric(p$np_diagnostics[[op_method]]$exact_zero_residual_count), numeric(1L))) else NA_real_,
      np_residual_abs_max = if (need_np) max(vapply(pieces, function(p)
        p$np_diagnostics[[op_method]]$residual_abs_max, numeric(1L))) else NA_real_,
      stringsAsFactors = FALSE)
  })
  parametric_residual <- numeric(n)
  for (piece in pieces) parametric_residual[piece$index] <- piece$observed$residual
  np_residual_diagnostics <- if (need_np) setNames(lapply(cfg$bandwidth_methods, function(method) {
    residual <- numeric(n)
    for (piece in pieces) residual[piece$index] <- piece$np_residual[[method]]
    residual_diagnostics(residual)
  }), cfg$bandwidth_methods) else NULL
  list(n = n, pieces = pieces, specs = specs,
       operator_methods = cfg$bandwidth_methods, summary = do.call(rbind, rows),
       parametric_residual = parametric_residual,
       parametric_residual_diagnostics = residual_diagnostics(parametric_residual),
       nonparametric_residual_diagnostics = np_residual_diagnostics)
}

fold_response_coefficient_delta <- function(piece, perturbation) {
  coefficient <- qr.coef(piece$observed$projection$qr, perturbation)
  matrix(as.numeric(coefficient), nrow = ncol(piece$observed$x),
         ncol = ncol(perturbation))
}

fold_operator_response_apply <- function(piece, method, perturbation) {
  if (piece$apply_backend == "blas")
    return(blas_dense_apply_cpp(piece$dense_operators[[method]], perturbation))
  fast_fold_operator_apply(piece$operators[[method]], perturbation,
                            backend = piece$apply_backend)
}

fold_response_statistic <- function(piece, component, smoothed_perturbation,
                                    coefficient_delta, cfg) {
  # xi* = X beta_hat + perturbation; beta* = beta_hat + coefficient_delta.
  # These are S xi* - X_eval beta* (canonical), or S xi* - S X beta*
  # (matched). Thus every response draw refits the fold null.
  evaluation_design <- piece$observed$evaluation_design
  if (component$null_fit == "canonical") {
    difference <- smoothed_perturbation - evaluation_design %*% coefficient_delta
    difference <- sweep(difference, 1L, component$smoothing_bias, "+")
  } else {
    smoothed_design <- evaluation_design + component$smoothing_bias_operator
    difference <- smoothed_perturbation - smoothed_design %*% coefficient_delta
  }
  piece$n * sqrt(piece$operators[[component$operator_method]]$h) *
    colSums(cfg$weighted_quadrature * difference^2)
}

fold_average_bootstrap_block <- function(state, omega, cfg,
                                         scheme = c("fold_parametric", "fold_nonparametric")) {
  scheme <- match.arg(scheme)
  answer <- matrix(0, ncol(omega), nrow(state$summary),
                   dimnames = list(NULL, state$summary$method))
  m <- length(cfg$quadrature_node)
  for (piece in state$pieces) {
    fold_omega <- omega[piece$index, , drop = FALSE]
    if (scheme == "fold_parametric") {
      perturbation <- piece$observed$residual * fold_omega
      delta <- fold_response_coefficient_delta(piece, perturbation)
      if (piece$apply_backend == "matrix") {
        smoothed <- as.matrix(piece$stacked %*% perturbation)
        for (j in seq_len(nrow(state$summary))) {
          component <- piece$observed$method[[state$summary$method[[j]]]]
          position <- match(component$operator_method, state$operator_methods)
          index <- seq.int((position - 1L) * m + 1L, position * m)
          answer[, j] <- answer[, j] + piece$weight * fold_response_statistic(
            piece, component, smoothed[index, , drop = FALSE], delta, cfg)
        }
      } else {
        for (op_method in state$operator_methods) {
          smoothed <- fold_operator_response_apply(piece, op_method, perturbation)
          for (j in which(state$summary$operator_method == op_method)) {
            component <- piece$observed$method[[state$summary$method[[j]]]]
            answer[, j] <- answer[, j] + piece$weight * fold_response_statistic(
              piece, component, smoothed, delta, cfg)
          }
        }
      }
    } else {
      # The canonical and matched comparisons share the same response draws.
      for (op_method in state$operator_methods) {
        perturbation <- piece$np_residual[[op_method]] * fold_omega
        delta <- fold_response_coefficient_delta(piece, perturbation)
        smoothed <- if (piece$apply_backend == "matrix")
          as.matrix(piece$operators[[op_method]]$S %*% perturbation) else
          fold_operator_response_apply(piece, op_method, perturbation)
        for (j in which(state$summary$operator_method == op_method)) {
          component <- piece$observed$method[[state$summary$method[[j]]]]
          answer[, j] <- answer[, j] + piece$weight * fold_response_statistic(
            piece, component, smoothed, delta, cfg)
        }
      }
    }
  }
  assert_true(all(is.finite(answer)) && all(answer >= 0),
              "Nonfinite or negative fold bootstrap statistic")
  answer
}

bootstrap_rows <- function(state, bootstrap, scheme, cfg) {
  max_exceedances <- floor(cfg$alpha * (cfg$bootstrap_draws + 1) - 1 + 1e-12)
  aligned_index <- cfg$bootstrap_draws - max_exceedances
  assert_true(aligned_index >= 1L && aligned_index <= cfg$bootstrap_draws,
              "Bootstrap size is incompatible with the requested alpha")
  rows <- lapply(seq_len(nrow(state$summary)), function(j) {
    result <- state$summary[j, , drop = FALSE]
    statistic <- result$statistic[[1L]]
    distribution <- bootstrap[, j]
    p_value <- (1 + sum(distribution >= statistic)) / (cfg$bootstrap_draws + 1)
    aligned_critical <- sort(distribution, partial = aligned_index)[[aligned_index]]
    reject <- p_value <= cfg$alpha
    assert_true(identical(reject, statistic > aligned_critical),
                "P-value and aligned critical-value decisions differ")
    result$bootstrap_scheme <- scheme
    result$residualization <- if (scheme == "fold_parametric") "parametric" else "nonparametric"
    result$multiplier_law <- "two_point"
    diagnostics <- if (scheme == "fold_parametric") state$parametric_residual_diagnostics else
      state$nonparametric_residual_diagnostics[[result$operator_method[[1L]]]]
    for (name in names(diagnostics)) result[[name]] <- unname(diagnostics[[name]])
    result$bootstrap_critical <- unname(quantile(distribution, 1 - cfg$alpha, type = 8))
    result$bootstrap_critical_aligned <- aligned_critical
    result$p_value <- p_value
    result$reject <- reject
    result
  })
  do.call(rbind, rows)
}

wild_bootstrap_tests <- function(y, a, degree, cfg, fold, fold_operators,
                                 fold_projections = NULL, random_blocks = NULL,
                                 multiplier_seed = NULL, np_fits = NULL,
                                 random_blocks_validated = FALSE) {
  assert_true(identical(cfg$multiplier_laws, "two_point"),
              "This experiment uses only two-point multipliers")
  assert_true(length(cfg$bootstrap_schemes) > 0L &&
                !anyDuplicated(cfg$bootstrap_schemes) &&
                all(cfg$bootstrap_schemes %in% c("fold_parametric", "fold_nonparametric")),
              "Unknown or repeated fold-bootstrap residualization")
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required for sparse bootstrap operators.", call. = FALSE)
  if (is.null(random_blocks)) {
    assert_true(!is.null(multiplier_seed), "A multiplier seed or cached blocks are required")
    random_blocks <- make_random_index_blocks(length(y), cfg, multiplier_seed)
  }
  if (!isTRUE(random_blocks_validated)) {
    validate_random_blocks(random_blocks, length(y), cfg)
  } else {
    # A caller reusing one common stream may perform the integer-entry scan
    # once. Keep inexpensive shape checks for each use of that validated stream.
    assert_true(length(random_blocks) > 0L &&
                  sum(vapply(random_blocks, ncol, integer(1L))) == cfg$bootstrap_draws &&
                  all(vapply(random_blocks, nrow, integer(1L)) == length(y)),
                "Previously validated random blocks have incompatible dimensions")
  }
  state <- fold_average_bootstrap_state(
    y, a, degree, cfg, fold, fold_operators, fold_projections, np_fits)
  bootstrap <- setNames(lapply(cfg$bootstrap_schemes, function(scheme)
    matrix(NA_real_, cfg$bootstrap_draws, nrow(state$summary),
           dimnames = list(NULL, state$summary$method))), cfg$bootstrap_schemes)
  start <- 1L
  for (index_block in random_blocks) {
    omega <- multiplier_from_index(index_block)
    stop_at <- start + ncol(index_block) - 1L
    for (scheme in cfg$bootstrap_schemes)
      bootstrap[[scheme]][start:stop_at, ] <- fold_average_bootstrap_block(
        state, omega, cfg, scheme)
    start <- stop_at + 1L
  }
  rows <- lapply(cfg$bootstrap_schemes, function(scheme)
    bootstrap_rows(state, bootstrap[[scheme]], scheme, cfg))
  residual <- state$parametric_residual
  list(result = do.call(rbind, rows), residual = residual,
       residual_diagnostics = state$parametric_residual_diagnostics,
       bootstrap = if (isTRUE(cfg$store_bootstrap_statistics)) bootstrap else NULL,
       state = if (isTRUE(cfg$store_bootstrap_statistics)) state else NULL)
}
