#!/usr/bin/env Rscript
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
source(file.path(root, "R", "utils.R"))
source(file.path(root, "config", "simulation_config.R"))
source(file.path(root, "R", "operators.R"))
source(file.path(root, "R", "wild_bootstrap.R"))
stopifnot(requireNamespace("Matrix", quietly = TRUE))

cfg <- simulation_config("test")
cfg$folds <- 3L
cfg$bootstrap_draws <- 39L
cfg$bootstrap_block <- 11L
cfg$bootstrap_schemes <- c("fold_parametric", "fold_nonparametric")
cfg$multiplier_laws <- "two_point"
cfg$bootstrap_apply_backend <- "matrix"
cfg$store_bootstrap_statistics <- TRUE
cfg$inference_min_neighbors <- 3L
sizes <- c(41L, 53L, 67L)
fold <- rep(seq_along(sizes), sizes)
a <- unlist(lapply(sizes, function(n) seq(-1.4, 1.4, length.out = n)))
y <- with_preserved_seed(561L, 2.5 + .3 * a + .2 * a^2 + .08 * a^3 + rnorm(length(a), sd = .35))
h <- setNames(c(.95, .8, .65, 1.03, .75, .31), cfg$bandwidth_methods)
fold_operators <- lapply(seq_len(cfg$folds), function(k)
  build_operator_set(a[fold == k], h, cfg, min_neighbors = 3L))
np_fits <- lapply(seq_len(cfg$folds), function(k) {
  ak <- a[fold == k]; yk <- y[fold == k]
  setNames(lapply(cfg$bandwidth_methods, function(method) {
    order <- if (method %in% cfg$pdr_methods) 2L else 1L
    smoother <- local_smoother_operator(ak, ak, h[[method]], order, cfg,
                                        min_neighbors = order + 1L)
    list(fitted = drop(smoother$S %*% yk), zero_count = 0L)
  }), cfg$bandwidth_methods)
})
blocks <- make_random_index_blocks(length(y), cfg, 12345L)
omega <- do.call(cbind, lapply(blocks, multiplier_from_index))
cfg_other_block <- cfg; cfg_other_block$bootstrap_block <- 17L
stopifnot(identical(do.call(cbind, blocks), do.call(cbind,
  make_random_index_blocks(length(y), cfg_other_block, 12345L))))
root5 <- sqrt(5); support <- (1 + c(-1, 1) * root5) / 2
probability <- c((5 + root5) / 10, (5 - root5) / 10)
stopifnot(abs(sum(support * probability)) < 1e-14,
          abs(sum(support^2 * probability) - 1) < 1e-14,
          abs(sum(support^3 * probability) - 1) < 1e-14)

max_error <- 0
for (degree in 0:2) {
  fit <- wild_bootstrap_tests(y, a, degree, cfg, fold, fold_operators,
                             random_blocks = blocks, np_fits = np_fits)
  specs <- inference_specs(degree, cfg)
  direct_observed <- numeric(nrow(specs))
  names(direct_observed) <- specs$method
  direct_bootstrap <- setNames(lapply(cfg$bootstrap_schemes, function(s)
    matrix(0, cfg$bootstrap_draws, nrow(specs), dimnames = list(NULL, specs$method))),
    cfg$bootstrap_schemes)
  for (k in seq_len(cfg$folds)) {
    index <- which(fold == k); nk <- length(index)
    X <- outer(a[index], 0:degree, `^`)
    Xeval <- outer(cfg$quadrature_node, 0:degree, `^`)
    qrx <- qr(X)
    beta <- qr.coef(qrx, y[index])
    null_fit <- drop(X %*% beta)
    for (j in seq_len(nrow(specs))) {
      method <- specs$method[[j]]
      op_method <- specs$operator_method[[j]]
      S <- fold_operators[[k]][[op_method]]$S
      comparator <- if (specs$null_fit[[j]] == "matched")
        drop(S %*% null_fit) else drop(Xeval %*% beta)
      difference <- drop(S %*% y[index]) - comparator
      direct_observed[[j]] <- direct_observed[[j]] +
        (nk / length(y)) * nk * sqrt(h[[op_method]]) *
        sum(cfg$weighted_quadrature * difference^2)
      for (scheme in cfg$bootstrap_schemes) {
        residual <- if (scheme == "fold_parametric") y[index] - null_fit else
          y[index] - np_fits[[k]][[op_method]]$fitted
        for (b in seq_len(cfg$bootstrap_draws)) {
          response <- null_fit + omega[index, b] * residual
          beta_star <- qr.coef(qrx, response)
          fitted_star <- if (specs$null_fit[[j]] == "matched")
            drop(S %*% drop(X %*% beta_star)) else drop(Xeval %*% beta_star)
          difference_star <- drop(S %*% response) - fitted_star
          direct_bootstrap[[scheme]][b, j] <- direct_bootstrap[[scheme]][b, j] +
            (nk / length(y)) * nk * sqrt(h[[op_method]]) *
            sum(cfg$weighted_quadrature * difference_star^2)
        }
      }
    }
  }
  for (scheme in cfg$bootstrap_schemes) {
    rows <- fit$result[fit$result$bootstrap_scheme == scheme, ]
    stopifnot(identical(rows$method, specs$method),
              max(abs(rows$statistic - direct_observed)) < 1e-10)
    error <- max(abs(fit$bootstrap[[scheme]] - direct_bootstrap[[scheme]]))
    max_error <- max(max_error, error)
    stopifnot(error < 1e-10)
    direct_p <- (1 + colSums(sweep(direct_bootstrap[[scheme]], 2L,
                                  direct_observed, `>=`))) / (cfg$bootstrap_draws + 1)
    stopifnot(identical(unname(rows$p_value), unname(direct_p)))
    for (j in seq_len(nrow(rows))) {
      residual <- numeric(length(y))
      for (k in seq_len(cfg$folds)) {
        index <- which(fold == k)
        X <- outer(a[index], 0:degree, `^`)
        residual[index] <- if (scheme == "fold_parametric")
          y[index] - drop(X %*% qr.coef(qr(X), y[index])) else
          y[index] - np_fits[[k]][[rows$operator_method[[j]]]]$fitted
      }
      expected <- residual_diagnostics(residual)
      stopifnot(max(abs(unlist(rows[j, names(expected)]) - expected)) < 1e-12)
    }
  }
  stopifnot(identical(fit$result$statistic[fit$result$bootstrap_scheme == "fold_parametric"],
                      fit$result$statistic[fit$result$bootstrap_scheme == "fold_nonparametric"]))
  # Each scheme must use its declared residual; the bootstrap distributions
  # differ even though the observed statistic and random-number indices agree.
  stopifnot(max(abs(fit$bootstrap$fold_parametric - fit$bootstrap$fold_nonparametric)) > 1e-5)
  # Unequal n_k verifies n_k/n weighting, not a hardcoded 1/K average.
  weights <- sizes / sum(sizes)
  unweighted <- colMeans(do.call(rbind, lapply(fit$state$pieces, function(piece)
    vapply(piece$observed$method, `[[`, numeric(1L), "statistic"))))
  stopifnot(any(abs(direct_observed - unweighted) > 1e-5),
            max(abs(weights - vapply(fit$state$pieces, `[[`, numeric(1L), "weight"))) < 1e-15)
  if (degree == 2L) {
    stopifnot(any(abs(direct_observed[c("ProxLL-raw", "ProxLL-US")] -
                      direct_observed[c("ProxLL-raw-matched", "ProxLL-US-matched")]) > 1e-5))
  }
  fit_prevalidated <- wild_bootstrap_tests(
    y, a, degree, cfg, fold, fold_operators, random_blocks = blocks,
    np_fits = np_fits, random_blocks_validated = TRUE)
  stopifnot(identical(fit$result, fit_prevalidated$result),
            identical(fit$bootstrap, fit_prevalidated$bootstrap))
}
diagnostic_fits <- np_fits
diagnostic_fits[[1L]][["PDR-PI1"]]$zero_index <- c(1L, 21L)
diagnostic_fits[[1L]][["PDR-PI1"]]$zero_count <- 2L
diagnostic_fits[[1L]][["PDR-PI1"]]$fitted[c(1L, 21L)] <- 0
diagnostic_fits[[1L]][["PDR-PI1"]]$max_condition <- 1e15
diagnostic_fits[[1L]][["PDR-PI1"]]$max_solve_condition <- 20
diagnostic_state <- fold_average_bootstrap_state(
  y, a, 2L, cfg, fold, fold_operators, np_fits = diagnostic_fits)
diagnostic_row <- diagnostic_state$summary[diagnostic_state$summary$method == "PDR-PI1", ]
stopifnot(diagnostic_row$np_zero_count == 2L,
          diagnostic_row$np_zero_inside_test_support == 1L,
          diagnostic_row$np_max_original_condition == 1e15,
          diagnostic_row$np_max_solve_condition == 20)
cat(sprintf("PASS: direct response/refit/integral equivalence, degrees 0:2, all methods, both residuals; max error %.3e\n", max_error))
cat("PASS: unequal fold weights, common two-point draws, multiplier moments, p-values and observed-statistic identity\n")
cat("PASS: method/residualization-specific residual diagnostics and prevalidated random-block path\n")
