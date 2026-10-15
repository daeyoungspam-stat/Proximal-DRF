#!/usr/bin/env Rscript
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(file_arg[[1L]])))
source(file.path(root, "R", "utils.R"))
source_simulation_files(root, compile_cpp = FALSE)
load_fold_fits(root)

cfg <- simulation_config("production")
cfg$folds <- 2L
cfg$quadrature_node <- c(-0.15, 0, 0.15)
cfg$fold_cv_filter <- TRUE
cfg$fold_cv_filter_rule <- "test-original-final-cloud-filter"
a <- rep(c(-1, -.8, -.6, .6, .8, 1), 2L)
fold <- rep(1:2, each = 6L)
hs <- c(.1, 1, 1.5)
risk <- matrix(c(0, 1, 2, 0, 2, 1), 3L, 2L,
               dimnames = list(NULL, c("first", "second")))
allowed <- list(hs, hs[c(1L, 3L)])
make_level <- function(risk) {
  chosen <- select_from_allowed(hs, risk, allowed)
  list(initial = hs, h0 = 1.5, cv1_ll = chosen,
       cv1_pdr = chosen, cv2_pdr = chosen,
       final = list(bandwidth = hs, risk_ll = risk, risk_pdr = risk))
}
level <- make_level(risk)
feas <- new_fold_feasibility(a, fold, cfg)
bad <- feas$check(.1, 2L)
stopifnot(!bad$feasible, length(bad$failures) == 2L,
          all(vapply(bad$failures, function(x)
            x$location == "quadrature" && x$local_n == 0L, logical(1L))))
for (h in hs[-1L]) for (degree in 1:2)
  stopifnot(feas$check(h, degree)$feasible)
stopifnot(feas$check(1 / log10(length(a)), 1L)$feasible)
# Exact rank deficiency at observed doses retains the original zero convention.
np <- nonparametric_fitted_values(a[fold == 1L], rep(0, 6L), .1, 2L, cfg)
stopifnot(np$zero_count == 6L, all(np$fitted == 0))
set.seed(97641)
rng_before <- .Random.seed
filtered <- filter_level_cv(level, a, fold, cfg, feas)
stopifnot(identical(.Random.seed, rng_before), identical(filtered$final, level$final),
          identical(filtered$initial, level$initial),
          identical(filtered$original_cv1_ll, level$cv1_ll),
          identical(filtered$admissibility$original$cv1_ll, level$cv1_ll))
for (selector in c("cv1_ll", "cv1_pdr", "cv2_pdr")) {
  stopifnot(identical(unname(filtered[[selector]]$h), c(1, 1.5)),
            identical(filtered[[selector]]$allowed, list(first = c(1, 1.5), second = 1.5)))
}
cat("PASS: sparse-fold candidates excluded; original risks and final candidate clouds preserved; NP rank-deficiency convention unchanged; no RNG consumption\n")

# If the original winner is feasible, the selected bandwidth and risk stay exact.
good_risk <- risk
good_risk[1L, ] <- 5
good_level <- make_level(good_risk)
good <- filter_level_cv(good_level, a, fold, cfg, feas)
for (selector in c("cv1_ll", "cv1_pdr", "cv2_pdr"))
  stopifnot(identical(good[[selector]]$h, good_level[[selector]]$h),
            identical(good[[selector]]$selected_risk,
                      good_level[[selector]]$selected_risk))
cat("PASS: feasible original winner and its CV risk unchanged\n")

# The LL admissible set uses BOTH the original h and h/log10(n), not only h.
pair_feas <- list(cache = new.env(parent = emptyenv()),
  check = function(h, degree) list(feasible = degree == 2L || h >= .98))
pair <- filter_level_cv(good_level, a, fold, cfg, pair_feas)
stopifnot(pair$cv1_ll$h[[1L]] == 1.5, pair$cv1_pdr$h[[1L]] == 1)
h_map <- bandwidth_map(1L, pair, list(list(h = .7)), length(a), cfg)
stopifnot(h_map[["ProxLL-US"]] == h_map[["ProxLL-raw"]] / log10(length(a)))
none <- list(cache = new.env(parent = emptyenv()),
             check = function(h, degree) list(feasible = FALSE))
empty_error <- tryCatch(filter_level_cv(level, a, fold, cfg, none), error = identity)
stopifnot(inherits(empty_error, "error"),
          grepl("No admissible CV candidate", conditionMessage(empty_error), fixed = TRUE))
cat("PASS: LL raw/US paired feasibility and exact ratio; empty admissible set errors explicitly\n")

# Unexpected implementation errors are never converted into inadmissible h.
probe_env <- new.env(parent = .GlobalEnv)
probe_env$local_smoother_operator <- function(...) stop("synthetic unexpected bug")
probe_factory <- new_fold_feasibility
environment(probe_factory) <- probe_env
unexpected <- tryCatch(probe_factory(a, fold, cfg)$check(1, 2L), error = identity)
stopifnot(inherits(unexpected, "error"),
          identical(conditionMessage(unexpected), "synthetic unexpected bug"))
cat("PASS: unexpected errors propagate without candidate masking\n")

# Inject inexpensive preparation components into the actual replication body.
# This verifies that PI1 still receives the ORIGINAL, unfiltered LL pilot even
# when the filtered LL estimator bandwidth changes; no bootstrap is executed.
runner_env <- new.env(parent = .GlobalEnv)
runner_env$captured_pilots <- numeric()
runner_env$validate_requested_design <- function(...) invisible(TRUE)
runner_env$simulation_conditions <- function(...) data.frame(condition = c("c1", "c2"))
runner_env$simulate_proximal_base <- function(n, seed, ...) list(A = a)
runner_env$outcome_matrix <- function(...) matrix(0, length(a), 2L)
runner_env$balanced_folds <- function(...) fold
runner_env$crossfit_all_conditions <- function(...) list(
  xi = matrix(0, length(a), 2L),
  key = data.frame(condition = c("c1", "c2"), scenario = c("first", "second")))
runner_env$select_level_bandwidths <- function(...) level
runner_env$level_plugin_bandwidth <- function(a, y, pilot_h, bounds, cfg) {
  runner_env$captured_pilots <- c(runner_env$captured_pilots, pilot_h)
  list(h = 1.5 + pilot_h, h_raw = 1.5 + pilot_h)
}
runner <- run_one_replication_impl
environment(runner) <- runner_env
old_prep <- Sys.getenv("PDR_PREP_ONLY", unset = NA_character_)
old_audit <- Sys.getenv("PDR_BANDWIDTH_AUDIT_DIR", unset = NA_character_)
Sys.setenv(PDR_PREP_ONLY = "1", PDR_BANDWIDTH_AUDIT_DIR = "")
audit <- tryCatch(runner(length(a), 373L, cfg, root, NULL), finally = {
  if (is.na(old_prep)) Sys.unsetenv("PDR_PREP_ONLY") else Sys.setenv(PDR_PREP_ONLY = old_prep)
  if (is.na(old_audit)) Sys.unsetenv("PDR_BANDWIDTH_AUDIT_DIR") else
    Sys.setenv(PDR_BANDWIDTH_AUDIT_DIR = old_audit)
})
stopifnot(identical(runner_env$captured_pilots, unname(level$cv1_ll$h)),
          all(vapply(seq_along(audit$original_h), function(j)
            identical(audit$original_h[[j]][c("PDR-PI1", "PDR-PI2")],
                      audit$filtered_h[[j]][c("PDR-PI1", "PDR-PI2")]), logical(1L))),
          audit$original_h[[1L]][["ProxLL-raw"]] !=
            audit$filtered_h[[1L]][["ProxLL-raw"]])
cat("PASS: actual preparation path preserves PI1 original pilot and both PI bandwidths despite LL reselection\n")

original_root <- file.path(dirname(root), "null_original")
original_env <- new.env(parent = .GlobalEnv)
sys.source(file.path(original_root, "R", "utils.R"), original_env)
sys.source(file.path(original_root, "config", "simulation_config.R"), original_env)
original_cfg <- original_env$simulation_config("production")
current_cfg <- simulation_config("production")
reps <- c(818L,819L,821L,835L,845L,889L,907L,910L,918L,925L,935L,983L,986L)
for (rep in reps) for (stream in c(1L, 2L, 500L))
  stopifnot(identical(seed_for(2500L, rep, stream, current_cfg$base_seed),
                     original_env$seed_for(2500L, rep, stream, original_cfg$base_seed)))
stopifnot(identical(deparse(body(seed_for)), deparse(body(original_env$seed_for))),
          identical(current_cfg$delta_grid, original_cfg$delta_grid),
          identical(current_cfg$bootstrap_draws, original_cfg$bootstrap_draws))
cat("PASS: all 39 failed-replication DGP/fold/multiplier seed mappings, departures and B unchanged from original ARC package\n")
