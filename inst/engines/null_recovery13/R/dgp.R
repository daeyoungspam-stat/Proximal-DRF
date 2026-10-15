dgp_parameters <- function() {
  list(
    sigma_x = 0.25,
    alpha0 = 0,
    alpha_x = c(0.25, 0.25),
    alpha_u = 0.15,
    sigma_a = sqrt(0.79),
    zeta0 = 0.25,
    zeta_t = 0.25,
    zeta_x = c(0.25, 0.25),
    zeta_u = 1,
    sigma_z = 0.5,
    mu0 = 0,
    mu_x = c(0.25, 0.25),
    mu_u = 1,
    sigma_w = 0.5,
    outcome_noise_halfwidth = 0.5,
    null_intercept = 2.5,
    null_linear = 0.4,
    null_quadratic = 0.25,
    outcome_proxy_scale = 0.4
  )
}

exposure_variance <- function(p = dgp_parameters()) {
  p$sigma_x^2 * sum(p$alpha_x^2) + p$alpha_u^2 + p$sigma_a^2
}

raw_alternative_direction <- function(a, degree) {
  degree <- as.integer(degree)
  if (!degree %in% 0:2) stop("Null degree must be 0, 1, or 2", call. = FALSE)
  variance_a <- exposure_variance()
  a^3 - 3 * variance_a * a
}

direction_normalizer <- function(degree, cfg) {
  raw <- raw_alternative_direction(cfg$quadrature_node, degree)
  value <- sqrt(weighted_integral(raw^2, cfg))
  assert_true(is.finite(value) && value > 0, "Invalid direction normalizer")
  value
}

alternative_direction <- function(a, degree, cfg) {
  raw_alternative_direction(a, degree) / direction_normalizer(degree, cfg)
}

null_curve <- function(a, degree, p = dgp_parameters()) {
  degree <- as.integer(degree)
  value <- rep(p$null_intercept, length(a))
  if (degree >= 1L) value <- value + p$null_linear * a
  if (degree >= 2L)
    value <- value + p$null_quadratic * (a^2 - exposure_variance(p))
  value
}

theta_truth <- function(a, degree, delta, cfg, p = dgp_parameters()) {
  null_curve(a, degree, p) + delta * alternative_direction(a, degree, cfg)
}

outcome_bridge_truth <- function(w, a, degree, delta, cfg,
                                 p = dgp_parameters()) {
  theta_truth(a, degree, delta, cfg, p) + p$outcome_proxy_scale * tanh(w)
}

simulate_proximal_base <- function(n, seed, keep_latent = FALSE) {
  p <- dgp_parameters()
  out <- with_preserved_seed(seed, {
    x1 <- rnorm(n, sd = p$sigma_x)
    x2 <- rnorm(n, sd = p$sigma_x)
    u <- rnorm(n)
    eps_a <- rnorm(n, sd = p$sigma_a)
    eps_z <- rnorm(n, sd = p$sigma_z)
    eps_w <- rnorm(n, sd = p$sigma_w)
    eps_y <- runif(n, -p$outcome_noise_halfwidth, p$outcome_noise_halfwidth)
    a <- p$alpha0 + p$alpha_x[[1L]] * x1 + p$alpha_x[[2L]] * x2 +
      p$alpha_u * u + eps_a
    z <- p$zeta0 + p$zeta_t * a + p$zeta_x[[1L]] * x1 +
      p$zeta_x[[2L]] * x2 + p$zeta_u * u + eps_z
    w_mean <- p$mu0 + p$mu_x[[1L]] * x1 + p$mu_x[[2L]] * x2 + p$mu_u * u
    w <- w_mean + eps_w
    data.frame(A = a, Z = z, W = w, X1 = x1, X2 = x2,
               U = u, W_mean = w_mean, eps_Y = eps_y,
               stringsAsFactors = FALSE)
  })
  if (!keep_latent) out[c("U", "W_mean")] <- NULL
  out
}

outcome_matrix <- function(base, conditions, cfg) {
  required <- c("condition", "degree", "delta")
  assert_true(all(required %in% names(conditions)), "Malformed simulation conditions")
  y <- vapply(seq_len(nrow(conditions)), function(j) {
    theta_truth(base$A, conditions$degree[[j]], conditions$delta[[j]], cfg) +
      dgp_parameters()$outcome_proxy_scale * tanh(base$W) + base$eps_Y
  }, numeric(nrow(base)))
  colnames(y) <- conditions$condition
  list(y = y, conditions = conditions)
}

treatment_bridge_moment_threshold <- function(p = dgp_parameters()) {
  variance_a <- exposure_variance(p)
  variance_residual <- p$sigma_z^2 +
    (p$zeta_u / p$alpha_u)^2 * p$sigma_a^2
  covariance_a_residual <- -(p$zeta_u / p$alpha_u) * p$sigma_a^2
  covariance <- matrix(c(variance_a, covariance_a_residual,
                         covariance_a_residual, variance_residual), 2L, 2L)
  denominator <- p$sigma_a^2 * p$zeta_u^2 + p$alpha_u^2 * p$sigma_z^2
  lambda <- p$alpha_u^2 / (2 * denominator)
  hessian <- diag(c(-1 / variance_a, 2 * lambda))
  positive_root <- max(Re(eigen(covariance %*% hessian,
                                only.values = TRUE)$values))
  1 / positive_root
}

treatment_bridge_truth <- function(z, a, x1, x2, p = dgp_parameters()) {
  variance_a <- exposure_variance(p)
  denominator <- p$sigma_a^2 * p$zeta_u^2 + p$alpha_u^2 * p$sigma_z^2
  lambda <- p$alpha_u^2 / (2 * denominator)
  integrability <- p$sigma_a^2 * p$zeta_u^2 / denominator
  mean_a_x <- p$alpha0 + p$alpha_x[[1L]] * x1 + p$alpha_x[[2L]] * x2
  residual <- z - p$zeta0 - p$zeta_t * a - p$zeta_x[[1L]] * x1 -
    p$zeta_x[[2L]] * x2 - (p$zeta_u / p$alpha_u) * (a - mean_a_x)
  log_q <- log(p$sigma_a / sqrt(variance_a)) + 0.5 * log(integrability) -
    a^2 / (2 * variance_a) + lambda * residual^2
  exp(log_q)
}

dgp_diagnostics <- function(base, cfg) {
  variance_a <- exposure_variance()
  matrix_data <- as.matrix(base[c("A", "Z", "W", "X1", "X2", "eps_Y")])
  orthogonality <- vapply(cfg$null_degrees, function(degree) {
    direction <- alternative_direction(base$A, degree, cfg)
    x <- outer(base$A, 0:degree, `^`)
    max(abs(colMeans(x * direction)))
  }, numeric(1L))
  c(
    min_a = min(base$A), max_a = max(base$A),
    mean_a = mean(base$A), variance_a = var(base$A),
    ks_normal_statistic = unname(suppressWarnings(
      ks.test(base$A, "pnorm", 0, sqrt(variance_a))$statistic)),
    fraction_test_support = mean(abs(base$A) <= max(abs(cfg$test_support))),
    expected_fraction_test_support =
      2 * pnorm(max(abs(cfg$test_support)) / sqrt(variance_a)) - 1,
    q_moment_threshold = treatment_bridge_moment_threshold(),
    max_abs_empirical_direction_moment = max(orthogonality),
    all_finite = as.numeric(all(is.finite(matrix_data)))
  )
}
