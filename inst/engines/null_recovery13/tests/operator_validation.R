#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(file_arg)) dirname(dirname(normalizePath(file_arg[[1L]]))) else getwd()
source(file.path(root, "R", "utils.R"))
source(file.path(root, "R", "operators.R"))
load_fold_fits(root)

cfg <- list(local_condition_limit = 1e12, local_reproduction_tolerance = 1e-8,
            inference_min_neighbors = 3L,
            bandwidth_methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                                  "ProxLL-raw", "ProxLL-US"),
            pdr_methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2"),
            methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                        "ProxLL-raw", "ProxLL-raw-matched",
                        "ProxLL-US", "ProxLL-US-matched"),
            bootstrap_schemes = c("fold_parametric", "fold_nonparametric"),
            multiplier_laws = "two_point")

reference_fit <- function(a, y, h, degree) {
  vapply(a, function(center) {
    u <- (a - center) / h
    active <- abs(u) < 1
    if (length(unique(a[active])) < degree + 1L) return(0)
    x <- outer(u[active], 0:degree, `^`)
    root_weight <- sqrt(kernel_epanechnikov(u[active]))
    fit <- lm.fit(x * root_weight, y[active] * root_weight)
    unname(fit$coefficients[[1L]])
  }, numeric(1L))
}
set.seed(93271)
a <- c(runif(180L, -1.2, 1.2), -4, 4)
y <- sin(a) + rnorm(length(a), sd = 0.2)
h <- 0.45
for (degree in 1:2) {
  fit <- nonparametric_fitted_values(a, y, h, degree, cfg)
  reference <- reference_fit(a, y, h, degree)
  stopifnot(max(abs(fit$fitted - reference)) < 1e-9,
            identical(fit$zero_index, c(181L, 182L)))
  permutation <- sample.int(length(a))
  shuffled <- nonparametric_fitted_values(a[permutation], y[permutation], h, degree, cfg)
  stopifnot(max(abs(shuffled$fitted[order(permutation)] - fit$fitted)) < 1e-11)
  for (power in 0:degree) {
    reproduced <- nonparametric_fitted_values(a, a^power, h, degree, cfg)
    active <- reproduced$unique_local_n >= degree + 1L
    stopifnot(max(abs(reproduced$fitted[active] - a[active]^power)) < 1e-10)
  }
  for (i in seq_along(a)[fit$unique_local_n >= degree + 1L]) {
    u <- (a - a[[i]]) / h
    x <- outer(u, 0:degree, `^`)
    gram <- crossprod(x, x * kernel_epanechnikov(u))
    stopifnot(abs(fit$condition[[i]] / kappa(gram, exact = TRUE) - 1) < 1e-8)
  }
}

# With b=h, the debiased local-linear intercept equals the local-quadratic
# intercept on full-rank windows, using the empirical bias-correction moment.
grid <- seq(-0.8, 0.8, length.out = 17L)
op <- local_smoother_operator(a, grid, h, 2L, cfg)
for (i in seq_along(grid)) {
  u <- (a - grid[[i]]) / h; weight <- kernel_epanechnikov(u)
  x1 <- cbind(1, u); x2 <- cbind(x1, u^2)
  d1 <- crossprod(x1, x1 * weight)
  beta1 <- solve(d1, crossprod(x1, y * weight))
  beta2 <- solve(crossprod(x2, x2 * weight), crossprod(x2, y * weight))
  c2 <- solve(d1, crossprod(x1, u^2 * weight))[[1L]]
  debiased <- beta1[[1L]] - c2 * beta2[[3L]]
  stopifnot(abs(debiased - beta2[[1L]]) < 1e-10,
            abs(drop(op$S[i, ] %*% y) - debiased) < 1e-10)
}
stopifnot(max(abs(as.matrix(op$S %*% cbind(1, a, a^2)) -
                    cbind(1, grid, grid^2))) < 1e-10)
for (degree in 1:2) {
  compact_op <- local_smoother_operator(a, grid, h, degree, cfg)
  compact <- t(vapply(seq_along(grid), function(i) {
    u <- (a - grid[[i]]) / h
    drop(outer(u, 0:(degree + 2L), `^`) %*%
           compact_op$polynomial_coefficients[i, ]) * (abs(u) < 1)
  }, numeric(length(a))))
  stopifnot(identical(compact_op$evaluation_grid, as.numeric(grid)),
            identical(compact_op$a, as.numeric(a)),
            max(abs(compact - as.matrix(compact_op$S))) < 1e-12)
}

# Rank is determined by distinct locations, not merely the observation count.
dup <- nonparametric_fitted_values(c(0, 0, 0.1, 0.1, 10), 1:5, 0.2, 2L, cfg)
stopifnot(identical(dup$zero_index, 1:5), all(dup$fitted == 0),
          identical(dup$unique_local_n, c(2L, 2L, 2L, 2L, 1L)))

# Pure basis scaling is not a failure and does not change the estimator.
close <- c(-1e-6, 0, 1e-6)
close_y <- 2 + 1e6 * close + 1e12 * close^2
scaled <- nonparametric_fitted_values(close, close_y, 10, 2L, cfg)
stopifnot(scaled$zero_count == 0L, scaled$max_condition > 1e12,
          scaled$max_solve_condition < 1e12,
          max(abs(scaled$fitted - close_y)) < 1e-8)

# A full-rank but numerically unresolved geometry is diagnosed, never zeroed.
bad <- tryCatch(nonparametric_fitted_values(c(0, 0.1, 0.1 + 1e-15), 1:3,
                                           0.2, 2L, cfg), error = identity)
stopifnot(inherits(bad, "pdr_sample_admissibility_error"))
stopifnot(nrow(bootstrap_scheme_specs(0L, "fold_parametric", cfg)) == 6L,
          nrow(bootstrap_scheme_specs(2L, "fold_nonparametric", cfg)) == 8L,
          bootstrap_cell_count(0L, cfg) == 12L,
          bootstrap_cell_count(2L, cfg) == 16L)

# Cache eviction changes storage only, not the smoother. Use a small test budget
# to exercise FIFO eviction and a direct uncached-versus-cached comparison.
cache <- new.env(parent = emptyenv())
cache$unrelated_metadata <- "preserve"
key1 <- operator_cache_key(2L, h)
key2 <- operator_cache_key(2L, h + .01)
key3 <- operator_cache_key(2L, h + .02)
entry_bytes <- as.numeric(object.size(op))
budget <- 2 * entry_bytes
operator_cache_put(cache, key1, op, budget)
operator_cache_put(cache, key2, op, budget)
stopifnot(identical(operator_cache_get(cache, key1, budget), op))
operator_cache_put(cache, key3, op, budget)
state <- operator_cache_metadata(cache)
stopifnot(!exists(key1, cache, inherits = FALSE),
          identical(state$order, c(key2, key3)), sum(state$bytes) <= budget,
          identical(cache$unrelated_metadata, "preserve"))
too_large_key <- operator_cache_key(1L, 1.7)
operator_cache_put(cache, too_large_key, op, entry_bytes - 1)
stopifnot(!exists(too_large_key, cache, inherits = FALSE),
          sum(operator_cache_metadata(cache)$bytes) <= entry_bytes - 1)
cfg$quadrature_node <- grid
h_methods <- setNames(rep(h, length(cfg$bandwidth_methods)), cfg$bandwidth_methods)
keys <- c(operator_cache_key(1L, h), operator_cache_key(2L, h))
uncached <- build_operator_set(a, h_methods, cfg)
fresh_cache <- new.env(parent = emptyenv())
cached1 <- build_operator_set(a, h_methods, cfg, fresh_cache, keys)
cached2 <- build_operator_set(a, h_methods, cfg, fresh_cache, keys)
stopifnot(identical(uncached, cached1), identical(uncached, cached2),
          sum(operator_cache_metadata(fresh_cache)$bytes) <= 32 * 1024^2)
cat("Operator validation passed: QR agreement, polynomial reproduction, b=h identity,\n",
    "rank-deficient zero convention, stable scaling, explicit numerical failure, method sets,\n",
    "bounded FIFO operator cache and cached-versus-uncached identity.\n", sep = "")
