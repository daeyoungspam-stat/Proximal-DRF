kernel_epanechnikov <- function(u) 0.75 * pmax(1 - u^2, 0)

stop_sample_admissibility <- function(message, ...) {
  condition <- c(list(message = message, call = NULL), list(...))
  class(condition) <- c("pdr_sample_admissibility_error", "error", "condition")
  stop(condition)
}

local_smoother_operator <- function(a, eval_grid, h, degree, cfg,
                                    min_neighbors = cfg$inference_min_neighbors) {
  n <- length(a); m <- length(eval_grid)
  assert_true(degree %in% 1:2, "Only local-linear and local-quadratic operators are used")
  assert_true(is.finite(h) && h > 0, "Invalid operator bandwidth")
  min_neighbors <- max(as.integer(min_neighbors), degree + 1L)
  smoother <- matrix(0, m, n)
  polynomial_coefficients <- matrix(0, m, degree + 3L)
  condition <- reproduction_error <- local_n <- numeric(m)
  powers <- 0:degree
  for (r in seq_len(m)) {
    u <- (a - eval_grid[[r]]) / h
    kh <- kernel_epanechnikov(u) / h
    active <- kh > 0
    if (sum(active) < min_neighbors)
      stop_sample_admissibility(
        "Too few observations in local test window",
        node_index = r, node = eval_grid[[r]], h = h, degree = degree,
        local_n = sum(active), min_neighbors = min_neighbors)
    w <- outer(u, powers, `^`)
    inverse <- local_inverse(crossprod(w, w * kh) / n, cfg$local_condition_limit)
    # On |u|<1 the exact smoother row is a polynomial of degree degree+2:
    # [.75/(n*h) e_0' D^{-1} r_degree(u)] (1-u^2).
    # The compact representation is used only for algebraically equivalent
    # response multiplication; the sparse operator and its checks are retained.
    coefficient <- .75 * inverse$inverse[, 1L] / (n * h)
    polynomial_coefficients[r, seq_len(degree + 1L)] <- coefficient
    polynomial_coefficients[r, seq.int(3L, degree + 3L)] <-
      polynomial_coefficients[r, seq.int(3L, degree + 3L)] - coefficient
    gamma <- drop(w %*% inverse$inverse[, 1L]) * kh
    smoother[r, ] <- gamma / n
    reproduction <- vapply(0:degree, function(power)
      sum(smoother[r, ] * (a - eval_grid[[r]])^power), numeric(1L))
    target <- c(1, rep(0, degree))
    condition[[r]] <- inverse$condition
    reproduction_error[[r]] <- max(abs(reproduction - target))
    local_n[[r]] <- sum(active)
  }
  if (max(reproduction_error) > cfg$local_reproduction_tolerance)
    stop(sprintf("Polynomial reproduction error %.3e exceeds tolerance",
                 max(reproduction_error)), call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required for sparse smoothing operators.", call. = FALSE)
  list(S = as(Matrix::Matrix(smoother, sparse = TRUE), "dgCMatrix"),
       h = h, degree = degree, condition = condition,
       reproduction_error = reproduction_error, local_n = local_n,
       evaluation_grid = as.numeric(eval_grid), a = as.numeric(a),
       polynomial_coefficients = polynomial_coefficients)
}

operator_cache_key <- function(degree, h) {
  paste0("d", as.integer(degree), "_h", sprintf("%a", as.numeric(h)))
}

operator_cache_metadata <- function(cache) {
  assert_true(is.environment(cache), "Operator cache must be an environment")
  state <- attr(cache, "pdr_operator_cache_fifo_v1", exact = TRUE)
  if (is.null(state)) {
    # Fresh simulation caches are empty. A prepopulated cache has no insertion
    # history; initialize its existing operator entries in deterministic order.
    keys <- ls(cache, all.names = TRUE)
    keys <- keys[grepl("^d[12]_h", keys)]
    keys <- keys[vapply(keys, function(key) {
      value <- get(key, envir = cache, inherits = FALSE)
      is.list(value) && !is.null(value$S) && !is.null(value$h) && !is.null(value$degree)
    }, logical(1L))]
    bytes <- setNames(vapply(keys, function(key)
      as.numeric(object.size(get(key, envir = cache, inherits = FALSE))), numeric(1L)), keys)
    state <- list(order = keys, bytes = bytes)
    attr(cache, "pdr_operator_cache_fifo_v1") <- state
  }
  state
}

operator_cache_trim <- function(cache, max_bytes = 32 * 1024^2) {
  assert_true(length(max_bytes) == 1L && is.finite(max_bytes) && max_bytes >= 0,
              "Invalid operator-cache byte limit")
  state <- operator_cache_metadata(cache)
  present <- vapply(state$order, exists, logical(1L), envir = cache, inherits = FALSE)
  state$order <- state$order[present]
  state$bytes <- state$bytes[state$order]
  while (length(state$order) && sum(state$bytes) > max_bytes) {
    oldest <- state$order[[1L]]
    rm(list = oldest, envir = cache)
    state$order <- state$order[-1L]
    state$bytes <- state$bytes[state$order]
  }
  attr(cache, "pdr_operator_cache_fifo_v1") <- state
  invisible(state)
}

operator_cache_get <- function(cache, key, max_bytes = 32 * 1024^2) {
  state <- operator_cache_trim(cache, max_bytes)
  if (key %in% state$order) get(key, envir = cache, inherits = FALSE) else NULL
}

operator_cache_put <- function(cache, key, value, max_bytes = 32 * 1024^2) {
  state <- operator_cache_trim(cache, max_bytes)
  if (key %in% state$order) {
    rm(list = key, envir = cache)
    state$order <- state$order[state$order != key]
    state$bytes <- state$bytes[state$order]
  }
  bytes <- as.numeric(object.size(value))
  if (bytes <= max_bytes) {
    while (length(state$order) && sum(state$bytes) + bytes > max_bytes) {
      oldest <- state$order[[1L]]
      rm(list = oldest, envir = cache)
      state$order <- state$order[-1L]
      state$bytes <- state$bytes[state$order]
    }
    assign(key, value, envir = cache)
    state$order <- c(state$order, key)
    state$bytes <- c(state$bytes, setNames(bytes, key))
  }
  attr(cache, "pdr_operator_cache_fifo_v1") <- state
  invisible(value)
}

build_operator_set <- function(a, h_by_method, cfg, cache = NULL,
                               cacheable_keys = NULL,
                               min_neighbors = cfg$inference_min_neighbors) {
  assert_true(identical(names(h_by_method), cfg$bandwidth_methods),
              "Bandwidth/method mismatch")
  answer <- lapply(cfg$bandwidth_methods, function(method) {
    degree <- if (method %in% cfg$pdr_methods) 2L else 1L
    key <- operator_cache_key(degree, h_by_method[[method]])
    use_cache <- !is.null(cache) && key %in% cacheable_keys
    if (use_cache) {
      cached <- operator_cache_get(cache, key)
      if (!is.null(cached)) return(cached)
    }
    value <- local_smoother_operator(
      a, cfg$quadrature_node, h_by_method[[method]], degree, cfg,
      min_neighbors = min_neighbors)
    if (use_cache) operator_cache_put(cache, key, value)
    value
  })
  names(answer) <- cfg$bandwidth_methods
  answer
}

null_design <- function(a, degree) outer(a, 0:as.integer(degree), `^`)

ols_projection <- function(x) {
  factor <- qr(x, LAPACK = FALSE)
  if (factor$rank != ncol(x)) stop("Rank-deficient parametric-null design", call. = FALSE)
  list(qr = factor, Q = qr.Q(factor, complete = FALSE))
}

null_projection <- function(a, degree) {
  x <- null_design(a, degree)
  c(list(x = x), ols_projection(x))
}

inference_specs <- function(degree, cfg) {
  degree <- as.integer(degree)
  answer <- data.frame(
    method = cfg$bandwidth_methods,
    operator_method = cfg$bandwidth_methods,
    null_fit = "canonical",
    stringsAsFactors = FALSE
  )
  if (degree == 2L) {
    answer <- rbind(
      answer,
      data.frame(
        method = c("ProxLL-raw-matched", "ProxLL-US-matched"),
        operator_method = c("ProxLL-raw", "ProxLL-US"),
        null_fit = "matched",
        stringsAsFactors = FALSE
      )
    )
  }
  assert_true(!anyDuplicated(answer$method), "Duplicate inference method label")
  assert_true(all(answer$method %in% cfg$methods), "Unknown inference method")
  assert_true(all(answer$operator_method %in% cfg$bandwidth_methods),
              "Unknown inference operator")
  answer[order(match(answer$method, cfg$methods)), , drop = FALSE]
}

active_methods <- function(degree, cfg) inference_specs(degree, cfg)$method

bootstrap_scheme_specs <- function(degree, scheme, cfg) {
  specs <- inference_specs(degree, cfg)
  if (scheme %in% c("fold_parametric", "fold_nonparametric")) return(specs)
  stop("Unknown bootstrap scheme: ", scheme, call. = FALSE)
}

load_fold_fits <- function(root = project_root()) {
  if (exists("epan_observation_fits_cpp", mode = "function", inherits = TRUE))
    return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE))
    stop("Package 'Rcpp' is required for observation-location local fits.",
         call. = FALSE)
  cache <- Sys.getenv("PDR_JOINT_CPP_CACHE", unset = "")
  if (!nzchar(cache)) cache <- file.path(root, ".cache", "Rcpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(root, "src", "fold_fits.cpp"), cacheDir = cache,
                  rebuild = FALSE, showOutput = FALSE)
  invisible(TRUE)
}

nonparametric_fitted_values <- function(a, y, h, degree, cfg) {
  assert_true(length(a) == length(y) && length(a) > 0L &&
                all(is.finite(a)) && all(is.finite(y)),
              "Invalid observation-location local-fit inputs")
  assert_true(degree %in% 1:2 && is.finite(h) && h > 0,
              "Invalid observation-location polynomial degree or bandwidth")
  if (!exists("epan_observation_fits_cpp", mode = "function", inherits = TRUE))
    load_fold_fits()
  result <- epan_observation_fits_cpp(
    as.numeric(a), as.numeric(y), h, as.integer(degree), cfg$local_condition_limit)
  if (length(result$bad_index)) {
    i <- result$bad_index[[1L]]
    stop_sample_admissibility(
      sprintf("Full-rank observation-location local design is numerically inadmissible (scaled kappa=%.3e)",
              result$solve_condition[[i]]),
      node_index = i, node = a[[i]], h = h, degree = degree,
      local_n = result$local_n[[i]], unique_local_n = result$unique_local_n[[i]],
      condition_number = result$condition[[i]],
      solve_condition_number = result$solve_condition[[i]], location = "observation")
  }
  result$bad_index <- NULL
  result$zero_count <- length(result$zero_index)
  full_rank <- result$unique_local_n >= degree + 1L
  result$max_condition <- if (any(full_rank))
    max(result$condition[full_rank]) else NA_real_
  result$max_solve_condition <- if (any(full_rank))
    max(result$solve_condition[full_rank]) else NA_real_
  result$min_local_n <- min(result$local_n)
  result
}

bootstrap_cell_count <- function(degree, cfg) {
  sum(vapply(cfg$bootstrap_schemes, function(scheme)
    nrow(bootstrap_scheme_specs(degree, scheme, cfg)), integer(1L))) *
    length(cfg$multiplier_laws)
}

expected_test_rows_per_replication <- function(cfg) {
  conditions <- simulation_conditions(cfg)
  sum(vapply(conditions$degree, function(degree)
    bootstrap_cell_count(degree, cfg), integer(1L))) * nrow(cfg$scenarios)
}

observed_test_components <- function(y, a, degree, operators, cfg,
                                     projection = NULL, specs = NULL) {
  if (is.null(projection)) projection <- null_projection(a, degree)
  if (is.null(specs)) specs <- inference_specs(degree, cfg)
  x <- projection$x
  evaluation_design <- null_design(cfg$quadrature_node, degree)
  beta <- qr.coef(projection$qr, y)
  residual <- drop(y - x %*% beta)
  fitted <- drop(x %*% beta)
  rows <- lapply(seq_len(nrow(specs)), function(j) {
    method <- specs$method[[j]]
    operator_method <- specs$operator_method[[j]]
    null_fit <- specs$null_fit[[j]]
    operator <- operators[[operator_method]]
    curve <- drop(operator$S %*% y)
    smoothed_design <- as.matrix(operator$S %*% x)
    smoothing_bias_operator <- smoothed_design - evaluation_design
    smoothing_bias <- drop(smoothing_bias_operator %*% beta)
    matched_difference <- drop(operator$S %*% residual)
    if (null_fit == "canonical") {
      comparator <- drop(evaluation_design %*% beta)
      difference <- curve - comparator
      represented <- matched_difference + smoothing_bias
    } else {
      comparator <- drop(operator$S %*% fitted)
      difference <- curve - comparator
      represented <- matched_difference
    }
    identity_error <- max(abs(difference - represented))
    statistic <- length(y) * sqrt(operator$h) *
      weighted_integral(represented^2, cfg)
    list(
      method = method, operator_method = operator_method, null_fit = null_fit,
      statistic = statistic, curve = curve, comparator = comparator,
      difference = represented, matched_difference = matched_difference,
      smoothing_bias_operator = smoothing_bias_operator,
      smoothing_bias = smoothing_bias, identity_error = identity_error
    )
  })
  names(rows) <- specs$method
  list(x = x, evaluation_design = evaluation_design, projection = projection,
       beta = beta, residual = residual, specs = specs, method = rows)
}
