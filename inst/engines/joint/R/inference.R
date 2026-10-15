local_polynomial_components <- function(a, y, eval, h, degree, derivative_order, cfg) {
  assert_true(degree >= derivative_order && derivative_order %in% 0:1,
              "Invalid local-polynomial derivative order")
  n <- length(a); u <- (a - eval) / h
  kh <- kernel_epanechnikov(u) / h
  active <- kh > 0
  if (sum(active) < cfg$inference_min_neighbors)
    stop("Too few observations in local inference window", call. = FALSE)
  powers <- 0:degree
  w <- outer(u, powers, `^`)
  design <- crossprod(w, w * kh) / n
  inverse <- local_inverse(design, cfg$local_condition_limit)
  moment <- colMeans(w * (kh * y))
  beta <- drop(inverse$inverse %*% moment)
  index <- derivative_order + 1L
  scale <- factorial(derivative_order) * h^(-derivative_order)
  gamma_weight <- scale * drop(w %*% inverse$inverse[, index]) * kh
  fitted <- drop(w %*% beta)
  design_correction <- gamma_weight * fitted
  estimate <- mean(gamma_weight * y)
  reproduction <- vapply(0:degree, function(power) {
    mean(gamma_weight * (a - eval)^power)
  }, numeric(1L))
  target <- numeric(degree + 1L)
  target[[derivative_order + 1L]] <- factorial(derivative_order)
  list(
    theta = estimate, Gamma = gamma_weight, gamma = design_correction,
    beta = beta, condition = inverse$condition, local_n = sum(active),
    reproduction = reproduction, reproduction_target = target,
    reproduction_error = max(abs(reproduction - target))
  )
}

method_polynomial_spec <- function(estimand, method) {
  pdr <- startsWith(method, "PDR-")
  if (estimand == "dose_response") {
    list(degree = if (pdr) 2L else 1L, derivative_order = 0L)
  } else if (estimand == "derivative") {
    # The literal derivative of ProxLL is the degree-one local-polynomial slope.
    list(degree = if (pdr) 3L else 1L, derivative_order = 1L)
  } else stop("Unknown estimand", call. = FALSE)
}

truth_for_estimand <- function(estimand, a) {
  if (estimand == "dose_response") theta_truth(a)
  else if (estimand == "derivative") theta_derivative_truth(a)
  else stop("Unknown estimand", call. = FALSE)
}

eif_at_grid <- function(a, xi, fold, h_delta, eval_grid, h,
                        degree, derivative_order, cfg) {
  n <- length(a); m <- length(eval_grid)
  estimate <- se <- phi_mean <- condition <- reproduction_error <- local_n <- numeric(m)
  phi <- matrix(NA_real_, n, m)
  for (j in seq_along(eval_grid)) {
    component <- local_polynomial_components(
      a, xi, eval_grid[[j]], h, degree, derivative_order, cfg
    )
    fold_mean_weight <- vapply(seq_len(cfg$folds), function(k)
      mean(component$Gamma[fold == k]), numeric(1L))
    product_measure_term <- fold_mean_weight[fold] * h_delta
    raw_phi <- component$Gamma * xi - component$gamma + product_measure_term
    phi_mean[[j]] <- mean(raw_phi)
    centered <- raw_phi - phi_mean[[j]]
    phi[, j] <- centered
    estimate[[j]] <- component$theta
    se[[j]] <- sqrt(sum(centered^2)) / n
    condition[[j]] <- component$condition
    reproduction_error[[j]] <- component$reproduction_error
    local_n[[j]] <- component$local_n
  }
  if (any(!is.finite(phi)) || any(!is.finite(se)) || any(se <= 0))
    stop("Non-finite or zero EIF variance", call. = FALSE)
  colnames(phi) <- format(eval_grid, trim = TRUE, scientific = FALSE)
  list(estimate = estimate, se = se, phi = phi, phi_mean = phi_mean,
       condition = condition, reproduction_error = reproduction_error,
       local_n = local_n)
}

evaluate_pointwise <- function(data, crossfit, scenario, method, estimand, h, cfg) {
  scenario_row <- cfg$scenarios[cfg$scenarios$scenario == scenario, , drop = FALSE]
  assert_true(nrow(scenario_row) == 1L, "Unknown nuisance scenario")
  outcome_model <- scenario_row$outcome_model[[1L]]
  spec <- method_polynomial_spec(estimand, method)
  result <- eif_at_grid(
    data$A, crossfit$xi[, scenario], crossfit$fold,
    crossfit$h_delta[, outcome_model], cfg$curve_grid, h,
    spec$degree, spec$derivative_order, cfg
  )
  z <- qnorm(1 - cfg$alpha / 2)
  curve <- data.frame(
    estimand = estimand, scenario = scenario, method = method,
    a = cfg$curve_grid,
    on_point_grid = grid_membership(cfg$curve_grid, cfg$point_grid),
    truth = truth_for_estimand(estimand, cfg$curve_grid),
    estimate = result$estimate, se = result$se,
    ci_lower = result$estimate - z * result$se,
    ci_upper = result$estimate + z * result$se,
    band_lower = NA_real_, band_upper = NA_real_, bandwidth = h,
    polynomial_degree = spec$degree, derivative_order = spec$derivative_order,
    phi_mean = result$phi_mean, max_condition = result$condition,
    reproduction_error = result$reproduction_error, local_n = result$local_n,
    stringsAsFactors = FALSE
  )
  list(curve = curve, phi = result$phi)
}

shared_gaussian_draws <- function(draws, dimension, seed) {
  with_preserved_seed(seed, matrix(rnorm(draws * dimension), nrow = draws))
}

gaussian_band_from_draws <- function(estimate, phi, alpha, z) {
  cross <- crossprod(phi)
  diagonal <- sqrt(diag(cross))
  if (any(!is.finite(diagonal)) || any(diagonal <= 0))
    stop("Gaussian band has a zero or non-finite EIF variance", call. = FALSE)
  raw <- cross / outer(diagonal, diagonal)
  raw <- (raw + t(raw)) / 2
  projection <- psd_correlation(raw)
  correlation <- projection$correlation
  eig <- eigen(correlation, symmetric = TRUE)
  root <- eig$vectors %*% (sqrt(pmax(eig$values, 0)) * t(eig$vectors))
  gaussian <- z %*% t(root)
  maxima <- apply(abs(gaussian), 1L, max)
  critical <- unname(quantile(maxima, 1 - alpha, type = 8))
  se <- diagonal / nrow(phi)
  list(
    critical = critical,
    lower = estimate - critical * se,
    upper = estimate + critical * se,
    min_eigen = projection$min_eigen,
    projected = projection$projected,
    projection_norm = projection$projection_norm
  )
}

attach_gaussian_band <- function(pointwise, cfg, z) {
  info <- gaussian_band_from_draws(
    pointwise$curve$estimate, pointwise$phi, cfg$alpha, z
  )
  curve <- pointwise$curve
  curve$band_lower <- info$lower; curve$band_upper <- info$upper
  list(
    curve = curve,
    band = data.frame(
      estimand = curve$estimand[[1L]], scenario = curve$scenario[[1L]],
      method = curve$method[[1L]], critical = info$critical,
      simultaneous_coverage = all(curve$truth >= curve$band_lower &
                                    curve$truth <= curve$band_upper),
      mean_width = mean(curve$band_upper - curve$band_lower),
      min_eigen = info$min_eigen, projected = info$projected,
      projection_norm = info$projection_norm,
      stringsAsFactors = FALSE
    )
  )
}
