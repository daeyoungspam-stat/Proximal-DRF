#!/usr/bin/env Rscript
# Portable numerical validation. Default: one core, no production-size timing.
# Optional: --cores 2; --benchmark adds the five-fold n=10000 block benchmark.
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
source(file.path(root, "R", "utils.R"))
source(file.path(root, "config", "simulation_config.R"))
source(file.path(root, "R", "operators.R"))
source(file.path(root, "R", "fast_apply.R"))
source(file.path(root, "R", "wild_bootstrap.R"))
args <- parse_cli()
cores <- cli_value(args, "cores", 1L, "integer")
benchmark <- cli_value(args, "benchmark", FALSE, "logical")
stopifnot(length(cores) == 1L, !is.na(cores), cores %in% 1:2)
if (.Platform$OS.type == "windows" && cores != 1L)
  stop("Use --cores 1 on Windows", call. = FALSE)
load_fast_fold_apply(root)
print(fast_fold_backend_diagnostics())

# The BLAS product must also respect zero dimensions and reject incompatible ones.
stopifnot(identical(dim(blas_dense_apply_cpp(matrix(numeric(), 0L, 3L), matrix(0, 3L, 2L))), c(0L, 2L)),
          identical(unname(blas_dense_apply_cpp(matrix(numeric(), 2L, 0L), matrix(numeric(), 0L, 3L))), matrix(0, 2L, 3L)),
          inherits(try(blas_dense_apply_cpp(matrix(0, 2L, 3L), matrix(0, 4L, 1L)), silent = TRUE), "try-error"))

jobs <- expand.grid(n = c(500L, 2000L), h = c(.02, .05, .2, .6, 2), degree = 1:2)
worker <- function(j) {
  cfg <- simulation_config("test")
  job <- jobs[j, ]
  # Stratification keeps every small-bandwidth local fit full rank.
  a <- seq(-1.1, 1.1, length.out = job$n)
  set.seed(1000L + j)
  response <- matrix(rnorm(job$n * 9L), job$n)
  response[, 1L] <- rep(c(-1, 1), length.out = job$n)
  response[, 2L] <- 1e8 * rep(c(-1, 1), length.out = job$n)
  response[, 3L] <- 1
  response[, 4L] <- a
  response[, 5L] <- a^2
  op <- local_smoother_operator(a, cfg$quadrature_node, job$h, job$degree,
                                cfg, min_neighbors = job$degree + 1L)
  ref <- as.matrix(op$S %*% response)
  got <- fast_fold_operator_apply(op, response, "blas")
  # This forward-error scale also handles deliberately cancelling 1e8 responses.
  scale <- pmax(1, max(Matrix::rowSums(abs(op$S))) * apply(abs(response), 2L, max))
  normalized_error <- max(sweep(abs(got - ref), 2L, scale, "/"))
  stopifnot(normalized_error < 2e-12)
  list(job = j, pid = Sys.getpid(), normalized_error = normalized_error)
}
products <- if (cores == 1L) lapply(seq_len(nrow(jobs)), worker) else
  parallel::mclapply(seq_len(nrow(jobs)), worker, mc.cores = cores,
                     mc.preschedule = TRUE, mc.set.seed = FALSE)
stopifnot(all(vapply(products, is.list, logical(1L))))
cat(sprintf("PASS: portable R-BLAS products (%d cores); maximum normalized error %.3e\n",
            cores, max(vapply(products, `[[`, numeric(1L), "normalized_error"))))

# Self-contained fixture: explicit responses, foldwise refits, unequal folds,
# all null degrees, all methods, and quadratic canonical/matched comparisons.
source(file.path(root, "tests", "bootstrap_equivalence.R"))
for (degree in 0:2) {
  cfg$bootstrap_apply_backend <- "matrix"
  reference <- wild_bootstrap_tests(y, a, degree, cfg, fold, fold_operators,
                                    random_blocks = blocks, np_fits = np_fits)
  cfg$bootstrap_apply_backend <- "blas"
  accelerated <- wild_bootstrap_tests(y, a, degree, cfg, fold, fold_operators,
                                      random_blocks = blocks, np_fits = np_fits)
  for (scheme in cfg$bootstrap_schemes)
    stopifnot(max(abs(reference$bootstrap[[scheme]] - accelerated$bootstrap[[scheme]])) < 1e-10)
  stopifnot(identical(reference$result$p_value, accelerated$result$p_value),
            identical(reference$result$reject, accelerated$result$reject),
            identical(reference$result$statistic, accelerated$result$statistic))
}
cat("PASS: R-BLAS bootstrap T, p-values and decisions match Matrix reference for all nulls and both residualizations\n")

if (isTRUE(benchmark)) {
  cfg <- simulation_config("production")
  cfg$bootstrap_draws <- 100L
  set.seed(9201L)
  fold <- rep(1:5, each = 2000L)
  a <- rnorm(10000L)
  y <- 2.5 + .4 * a + .25 * a^2 + rnorm(10000L, sd = .35)
  h <- setNames(c(.2, .24, .35, .5, .35, .0875), cfg$bandwidth_methods)
  fold_operators <- lapply(1:5, function(k)
    build_operator_set(a[fold == k], h, cfg, min_neighbors = 3L))
  # Fixed arbitrary response vectors compare algorithms; no simulation DGP run.
  np_fits <- lapply(1:5, function(k) setNames(lapply(cfg$bandwidth_methods, function(method)
    list(fitted = 2.5 + .4 * a[fold == k] + .25 * a[fold == k]^2,
         zero_index = integer(), zero_count = 0L)), cfg$bandwidth_methods))
  omega <- multiplier_from_index(make_random_index_blocks(10000L, cfg, 8674L)[[1L]])
  timing <- setNames(numeric(2L), c("matrix", "blas"))
  values <- list()
  for (backend in names(timing)) {
    cfg$bootstrap_apply_backend <- backend
    state <- fold_average_bootstrap_state(y, a, 2L, cfg, fold, fold_operators, np_fits = np_fits)
    values[[backend]] <- lapply(cfg$bootstrap_schemes, function(scheme)
      fold_average_bootstrap_block(state, omega, cfg, scheme))
    timing[[backend]] <- system.time(for (iteration in 1:3)
      for (scheme in cfg$bootstrap_schemes)
        fold_average_bootstrap_block(state, omega, cfg, scheme))[["elapsed"]] / 3
  }
  stopifnot(max(abs(unlist(values$matrix) - unlist(values$blas))) < 1e-10)
  print(timing)
  cat(sprintf("Production-size block time ratio Matrix / R-BLAS: %.2f\n",
              timing[["matrix"]] / timing[["blas"]]))
}
