# Diagnostic rule: original final candidate clouds and CV risks are retained.
# Only design-based feasibility is added; no response or bootstrap result is
# used by the feasibility predicate. PI1's original pilot is preserved.
new_fold_feasibility <- function(a, fold, cfg) {
  stopifnot(length(a) == length(fold), all(is.finite(a)))
  indices <- lapply(seq_len(cfg$folds), function(k) which(fold == k))
  cache <- new.env(parent = emptyenv())
  check <- function(h, degree) {
    key <- operator_cache_key(degree, h)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    failures <- list()
    for (k in seq_along(indices)) {
      ak <- a[indices[[k]]]
      for (location in c("quadrature", "observation")) {
        result <- tryCatch({
          if (location == "quadrature") {
            local_smoother_operator(ak, cfg$quadrature_node, h, degree, cfg,
                                   min_neighbors = cfg$fold_inference_min_neighbors)
          } else {
            # Preserve the original zero-fit convention on rank-deficient
            # observed-dose windows; only original raised errors are excluded.
            nonparametric_fitted_values(ak, numeric(length(ak)), h, degree, cfg)
          }
          NULL
        }, error = function(e) e)
        if (inherits(result, "error")) {
          known <- inherits(result, "pdr_sample_admissibility_error") ||
            grepl("^Ill-conditioned local design|^Polynomial reproduction error|singular",
                  conditionMessage(result), ignore.case = TRUE)
          if (!known) stop(result)
          failures[[length(failures) + 1L]] <- c(list(fold = k,
            location = location, error = conditionMessage(result)),
            as.list(result[setdiff(names(result), c("message", "call", "location"))]))
        }
      }
    }
    answer <- list(h = h, degree = degree, feasible = !length(failures),
                   failures = failures)
    assign(key, answer, cache)
    answer
  }
  list(check = check, cache = cache)
}

filter_level_cv <- function(level, a, fold, cfg, feasibility = NULL) {
  if (is.null(feasibility)) feasibility <- new_fold_feasibility(a, fold, cfg)
  original <- level
  h <- level$final$bandwidth
  pdr_needed <- unique(c(unlist(level$cv1_pdr$allowed),
                        unlist(level$cv2_pdr$allowed)))
  ll_needed <- unique(unlist(level$cv1_ll$allowed))
  matches <- function(value, cloud) any(abs(log(value / cloud)) < 1e-12)
  pdr_ok <- ll_ok <- rep(FALSE, length(h))
  for (j in seq_along(h)) {
    if (matches(h[[j]], pdr_needed))
      pdr_ok[[j]] <- feasibility$check(h[[j]], 2L)$feasible
    if (matches(h[[j]], ll_needed)) {
      raw_ok <- feasibility$check(h[[j]], 1L)$feasible
      us_ok <- feasibility$check(h[[j]] / log10(length(a)), 1L)$feasible
      ll_ok[[j]] <- raw_ok && us_ok
    }
  }
  select <- function(risk, allowed, ok) {
    risk[!ok, ] <- Inf
    allowed <- lapply(allowed, function(values)
      values[vapply(values, function(v) any(ok & abs(log(h / v)) < 1e-12), logical(1))])
    select_from_allowed(h, risk, allowed)
  }
  level$cv1_ll <- select(level$final$risk_ll, level$cv1_ll$allowed, ll_ok)
  level$cv1_pdr <- select(level$final$risk_pdr, level$cv1_pdr$allowed, pdr_ok)
  level$cv2_pdr <- select(level$final$risk_pdr, level$cv2_pdr$allowed, pdr_ok)
  level$original_cv1_ll <- original$cv1_ll
  level$admissibility <- list(rule = cfg$fold_cv_filter_rule,
    original = original[c("cv1_ll", "cv1_pdr", "cv2_pdr")],
    mask = data.frame(h = h, pdr_feasible = pdr_ok, ll_pair_feasible = ll_ok),
    checks = as.list(feasibility$cache))
  level
}
