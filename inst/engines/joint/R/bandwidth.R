load_fast_cv <- function(root = project_root()) {
  needed <- c("epan_level_loocv_cpp", "epan_loo_density_cpp",
              "epan_derivative_loocv_cpp")
  if (all(vapply(needed, exists, logical(1), mode = "function", inherits = TRUE)))
    return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE))
    stop("Package 'Rcpp' is required for accelerated cross-validation.", call. = FALSE)
  cache <- Sys.getenv("PDR_JOINT_CPP_CACHE", unset = "")
  if (!nzchar(cache)) cache <- file.path(root, ".cache", "Rcpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(root, "src", "fast_cv.cpp"), cacheDir = cache,
                  rebuild = FALSE, showOutput = FALSE)
  invisible(TRUE)
}

robust_exposure_scale <- function(a) {
  scale <- min(sd(a), IQR(a) / 1.349)
  if (!is.finite(scale) || scale <= 0) scale <- sd(a)
  assert_true(is.finite(scale) && scale > 0, "Invalid exposure scale")
  scale
}

reference_bandwidth <- function(a, exponent, multiplier = 1.06) {
  multiplier * robust_exposure_scale(a) * length(a)^exponent
}

pi2_bandwidth <- function(a, exponent, cfg) {
  h <- cfg$pi2_multiplier * sd(a) * length(a)^exponent
  assert_true(is.finite(h) && h > 0, "Invalid PI2 bandwidth")
  unname(h)
}

cv2_candidates <- function(center, cfg) {
  h <- center * 2^(cfg$cv2_j * cfg$cv2_log2_step)
  assert_true(all(is.finite(h)) && all(h > 0) && all(diff(h) > 0),
              "Invalid CV2 grid")
  h
}

level_loocv <- function(a, xi, bandwidths, cfg) {
  trim <- quantile(a, c(cfg$level_cv_trim, 1 - cfg$level_cv_trim), names = FALSE)
  epan_level_loocv_cpp(
    as.numeric(a), as.matrix(xi), sort(unique(as.numeric(bandwidths))),
    trim[[1L]], trim[[2L]], cfg$cv_min_neighbors, cfg$cv_leverage_limit
  )
}

derivative_density <- function(a, cfg) {
  g <- reference_bandwidth(a, -1 / 5)
  raw <- epan_loo_density_cpp(as.numeric(a), g)
  floor_hit <- !is.finite(raw) | raw < cfg$derivative_density_floor
  used <- pmax(raw, cfg$derivative_density_floor)
  list(raw = raw, used = used, inverse = 1 / used, bandwidth = g,
       floor_hits = sum(floor_hit), floor_fraction = mean(floor_hit),
       min_raw = min(raw), median_raw = median(raw))
}

derivative_loocv <- function(a, xi, inverse_density, bandwidths, cfg) {
  epan_derivative_loocv_cpp(
    as.numeric(a), as.matrix(xi), as.numeric(inverse_density),
    sort(unique(as.numeric(bandwidths))),
    cfg$derivative_cv_domain[[1L]], cfg$derivative_cv_domain[[2L]],
    cfg$derivative_cv_flat_domain[[1L]], cfg$derivative_cv_flat_domain[[2L]],
    cfg$cv_min_neighbors
  )
}

adjacent_refinement <- function(initial, risk) {
  initial <- sort(unique(initial)); risk <- as.matrix(risk)
  lapply(seq_len(ncol(risk)), function(column) {
    finite <- which(is.finite(risk[, column]))
    if (!length(finite)) return(initial)
    selected <- finite[[which.min(risk[finite, column])]]
    addition <- numeric()
    if (selected > 1L)
      addition <- c(addition, sqrt(initial[[selected - 1L]] * initial[[selected]]))
    if (selected < length(initial))
      addition <- c(addition, sqrt(initial[[selected]] * initial[[selected + 1L]]))
    sort(unique(c(initial, addition)))
  })
}

select_from_allowed <- function(bandwidth, risk, allowed) {
  risk <- as.matrix(risk)
  scenarios <- colnames(risk)
  assert_true(length(allowed) == length(scenarios),
              "CV admissible-set/scenario dimension mismatch")
  names(allowed) <- scenarios
  selected <- selected_risk <- second_risk <- gap <- normalized_gap <- numeric(ncol(risk))
  selected_index <- integer(ncol(risk)); lower_hit <- upper_hit <- logical(ncol(risk))
  for (column in seq_len(ncol(risk))) {
    mask <- vapply(bandwidth, function(h)
      any(abs(log(h / allowed[[column]])) < 1e-12), logical(1))
    finite <- which(mask & is.finite(risk[, column]))
    if (!length(finite)) stop("No admissible CV candidate for ", scenarios[[column]], call. = FALSE)
    order <- finite[order(risk[finite, column], bandwidth[finite])]
    index <- order[[1L]]
    selected[[column]] <- bandwidth[[index]]
    selected_index[[column]] <- index
    selected_risk[[column]] <- risk[index, column]
    second_risk[[column]] <- if (length(order) >= 2L) risk[order[[2L]], column] else NA_real_
    gap[[column]] <- second_risk[[column]] - selected_risk[[column]]
    normalized_gap[[column]] <- gap[[column]] / (abs(selected_risk[[column]]) + 1e-12)
    # Boundary diagnostics refer to the admissible finite-risk set.  A formally
    # listed candidate with a failed local fit cannot be a selectable endpoint.
    allowed_range <- range(bandwidth[finite])
    lower_hit[[column]] <- abs(log(selected[[column]] / allowed_range[[1L]])) < 1e-12
    upper_hit[[column]] <- abs(log(selected[[column]] / allowed_range[[2L]])) < 1e-12
  }
  names(selected) <- names(selected_index) <- names(selected_risk) <- names(second_risk) <-
    names(gap) <- names(normalized_gap) <- names(lower_hit) <- names(upper_hit) <- scenarios
  list(h = selected, selected_index = selected_index, selected_risk = selected_risk,
       second_risk = second_risk, risk_gap = gap, normalized_risk_gap = normalized_gap,
       lower_hit = lower_hit, upper_hit = upper_hit,
       boundary_hit = lower_hit | upper_hit, allowed = allowed)
}

select_level_bandwidths <- function(a, xi, cfg) {
  base <- reference_bandwidth(a, -1 / 5)
  initial <- base * cfg$cv1_multipliers
  first <- level_loocv(a, xi, initial, cfg)
  ll_allowed <- if (isTRUE(cfg$cv1_refine))
    adjacent_refinement(first$bandwidth, first$risk_ll) else rep(list(initial), ncol(xi))
  pdr_allowed <- if (isTRUE(cfg$cv1_refine))
    adjacent_refinement(first$bandwidth, first$risk_pdr) else rep(list(initial), ncol(xi))
  h0 <- pi2_bandwidth(a, -1 / 5, cfg)
  cv2 <- cv2_candidates(h0, cfg)
  final_candidates <- sort(unique(c(unlist(ll_allowed), unlist(pdr_allowed), cv2)))
  final <- level_loocv(a, xi, final_candidates, cfg)
  cv1_ll <- select_from_allowed(final$bandwidth, final$risk_ll, ll_allowed)
  cv1_pdr <- select_from_allowed(final$bandwidth, final$risk_pdr, pdr_allowed)
  cv2_allowed <- rep(list(cv2), ncol(xi))
  cv2_pdr <- select_from_allowed(final$bandwidth, final$risk_pdr, cv2_allowed)
  list(base = base, initial = initial, h0 = h0, cv2_candidates = cv2,
       cv1_ll = cv1_ll, cv1_pdr = cv1_pdr, cv2_pdr = cv2_pdr,
       final = final)
}

select_derivative_bandwidths <- function(a, xi, cfg) {
  density <- derivative_density(a, cfg)
  base <- reference_bandwidth(a, cfg$derivative_rate_exponent)
  initial <- base * cfg$cv1_multipliers
  first <- derivative_loocv(a, xi, density$inverse, initial, cfg)
  ll_allowed <- if (isTRUE(cfg$cv1_refine))
    adjacent_refinement(first$bandwidth, first$risk_ll) else rep(list(initial), ncol(xi))
  pdr_allowed <- if (isTRUE(cfg$cv1_refine))
    adjacent_refinement(first$bandwidth, first$risk_pdr) else rep(list(initial), ncol(xi))
  h0 <- pi2_bandwidth(a, cfg$derivative_rate_exponent, cfg)
  cv2 <- cv2_candidates(h0, cfg)
  final_candidates <- sort(unique(c(unlist(ll_allowed), unlist(pdr_allowed), cv2)))
  final <- derivative_loocv(a, xi, density$inverse, final_candidates, cfg)
  cv1_ll <- select_from_allowed(final$bandwidth, final$risk_ll, ll_allowed)
  cv1_pdr <- select_from_allowed(final$bandwidth, final$risk_pdr, pdr_allowed)
  cv2_pdr <- select_from_allowed(final$bandwidth, final$risk_pdr, rep(list(cv2), ncol(xi)))
  list(base = base, initial = initial, h0 = h0, cv2_candidates = cv2,
       density = density, cv1_ll = cv1_ll, cv1_pdr = cv1_pdr,
       cv2_pdr = cv2_pdr, final = final)
}

rolling_conditional_variance <- function(a, y, power = 0.5) {
  n <- length(a); order <- order(a); sorted_y <- y[order]
  k <- max(15L, as.integer(ceiling(n^power)))
  if (k %% 2L == 0L) k <- k + 1L
  k <- min(k, if (n %% 2L == 1L) n else n - 1L)
  half <- k %/% 2L
  cumulative <- c(0, cumsum(sorted_y)); cumulative2 <- c(0, cumsum(sorted_y^2))
  sorted_answer <- numeric(n)
  for (i in seq_len(n)) {
    left <- max(1L, i - half); right <- min(n, left + k - 1L)
    left <- max(1L, right - k + 1L); count <- right - left + 1L
    total <- cumulative[[right + 1L]] - cumulative[[left]]
    total2 <- cumulative2[[right + 1L]] - cumulative2[[left]]
    sorted_answer[[i]] <- max((total2 - total^2 / count) / max(1, count - 1L), 1e-10)
  }
  answer <- numeric(n); answer[order] <- sorted_answer; answer
}

kernel_epanechnikov <- function(u) 0.75 * pmax(1 - u^2, 0)

level_pi_components <- function(a, y, eval, pilot_h, sigma2, cfg) {
  u <- (a - eval) / pilot_h; kh <- kernel_epanechnikov(u) / pilot_h
  active <- kh > 0
  if (sum(active) < cfg$cv_min_neighbors) stop("Too few level PI observations")
  w1 <- cbind(1, u); w2 <- cbind(1, u, u^2); n <- length(a)
  d1 <- crossprod(w1, w1 * kh) / n; d2 <- crossprod(w2, w2 * kh) / n
  i1 <- local_inverse(d1, cfg$local_condition_limit)
  i2 <- local_inverse(d2, cfg$local_condition_limit)
  m2 <- colMeans(w2 * (kh * y)); beta2 <- drop(i2$inverse %*% m2)
  theta2 <- 2 * beta2[[3L]] / pilot_h^2
  bias_vector <- colMeans(w1 * (kh * (0.5 * theta2 * u^2)))
  bias <- drop(i1$inverse[1L, , drop = FALSE] %*% bias_vector)
  middle <- crossprod(w1, w1 * (kh^2 * sigma2)) / n
  raw_variance <- drop(i1$inverse[1L, , drop = FALSE] %*% middle %*%
                         i1$inverse[, 1L, drop = FALSE])
  variance <- pilot_h * raw_variance
  c(B = bias, V = variance, theta2 = theta2,
    condition_primary = i1$condition, condition_pilot = i2$condition,
    neighbors = sum(active))
}

level_plugin_bandwidth <- function(a, y, h_ll_cv, bounds, cfg) {
  pilot_h <- min(max(cfg$level_pi_pilot_multiplier * h_ll_cv, bounds[[1L]]), bounds[[2L]])
  sigma2 <- rolling_conditional_variance(a, y, cfg$variance_knn_power)
  components <- t(vapply(cfg$curve_grid, function(eval) {
    tryCatch(level_pi_components(a, y, eval, pilot_h, sigma2, cfg),
             error = function(e) rep(NA_real_, 6L))
  }, numeric(6L)))
  colnames(components) <- c("B", "V", "theta2", "condition_primary",
                            "condition_pilot", "neighbors")
  valid <- apply(components, 1L, function(value) all(is.finite(value))) & components[, "V"] > 0
  if (sum(valid) < max(5L, length(cfg$curve_grid) %/% 2L))
    stop("Level PI pilot failed on too many grid points", call. = FALSE)
  ib <- trapz_grid(cfg$curve_grid[valid], components[valid, "B"]^2)
  iv <- trapz_grid(cfg$curve_grid[valid], components[valid, "V"])
  if (!is.finite(ib) || ib <= 1e-12 || !is.finite(iv) || iv <= 0)
    stop("Invalid level PI AMISE components", call. = FALSE)
  raw <- length(a)^(-1 / 5) * (iv / (4 * ib))^(1 / 5)
  used <- min(max(raw, bounds[[1L]]), bounds[[2L]])
  list(h = used, h_raw = raw, pilot_h = pilot_h, clipped = abs(log(raw / used)) > 1e-12,
       integrated_bias2 = ib, integrated_variance = iv, components = components)
}

derivative_pi_components <- function(a, y, eval, pilot_h, sigma2, cfg) {
  u <- (a - eval) / pilot_h; kh <- kernel_epanechnikov(u) / pilot_h
  active <- kh > 0
  if (sum(active) < cfg$cv_min_neighbors) stop("Too few derivative PI observations")
  w2 <- cbind(1, u, u^2); w3 <- cbind(1, u, u^2, u^3); n <- length(a)
  d2 <- crossprod(w2, w2 * kh) / n; d3 <- crossprod(w3, w3 * kh) / n
  i2 <- local_inverse(d2, cfg$local_condition_limit)
  i3 <- local_inverse(d3, cfg$local_condition_limit)
  m3 <- colMeans(w3 * (kh * y)); beta3 <- drop(i3$inverse %*% m3)
  theta3 <- 6 * beta3[[4L]] / pilot_h^3
  u3 <- colMeans(w2 * (kh * u^3))
  c3 <- drop(i2$inverse[2L, , drop = FALSE] %*% u3)
  bias <- c3 * theta3 / 6
  middle <- crossprod(w2, w2 * (kh^2 * sigma2)) / n
  raw_variance <- drop(i2$inverse[2L, , drop = FALSE] %*% middle %*%
                         i2$inverse[, 2L, drop = FALSE])
  variance <- pilot_h * raw_variance
  c(B = bias, V = variance, theta3 = theta3, c3 = c3,
    condition_primary = i2$condition, condition_pilot = i3$condition,
    neighbors = sum(active))
}

derivative_plugin_bandwidth <- function(a, y, bounds, cfg) {
  pilot_h <- cfg$derivative_pi_pilot_multiplier * robust_exposure_scale(a) *
    length(a)^cfg$derivative_pi_pilot_exponent
  sigma2 <- rolling_conditional_variance(a, y, cfg$variance_knn_power)
  components <- t(vapply(cfg$curve_grid, function(eval) {
    tryCatch(derivative_pi_components(a, y, eval, pilot_h, sigma2, cfg),
             error = function(e) rep(NA_real_, 7L))
  }, numeric(7L)))
  colnames(components) <- c("B", "V", "theta3", "c3", "condition_primary",
                            "condition_pilot", "neighbors")
  valid <- apply(components, 1L, function(value) all(is.finite(value))) & components[, "V"] > 0
  if (sum(valid) < max(5L, length(cfg$curve_grid) %/% 2L))
    stop("Derivative PI pilot failed on too many grid points", call. = FALSE)
  ib <- trapz_grid(cfg$curve_grid[valid], components[valid, "B"]^2)
  iv <- trapz_grid(cfg$curve_grid[valid], components[valid, "V"])
  if (!is.finite(ib) || ib <= 1e-12 || !is.finite(iv) || iv <= 0)
    stop("Invalid derivative PI AMISE components", call. = FALSE)
  raw <- (3 * iv / (4 * length(a) * ib))^(1 / 7)
  used <- min(max(raw, bounds[[1L]]), bounds[[2L]])
  list(h = used, h_raw = raw, pilot_h = pilot_h, clipped = abs(log(raw / used)) > 1e-12,
       integrated_bias2 = ib, integrated_variance = iv, components = components)
}
