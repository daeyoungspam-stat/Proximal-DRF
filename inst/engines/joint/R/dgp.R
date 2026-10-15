dgp_parameters <- function() {
  list(
    sigma_x = 0.25,
    alpha0 = 0,
    alpha_x = c(0.25, 0.25),
    alpha_u = 0.15,
    sigma_a = sqrt(0.79),
    zeta0 = 0.25,
    zeta_a = 0.25,
    zeta_x = c(0.25, 0.25),
    zeta_u = 1,
    sigma_z = 0.5,
    mu0 = 0,
    mu_x = c(0.25, 0.25),
    mu_u = 1,
    sigma_w = 0.5,
    dose_shift = 0.2,
    outcome_noise_halfwidth = 0.5
  )
}

theta_truth <- function(a, p = dgp_parameters()) 2.5 + sin(a + p$dose_shift)
theta_derivative_truth <- function(a, p = dgp_parameters()) cos(a + p$dose_shift)
theta_second_truth <- function(a, p = dgp_parameters()) -sin(a + p$dose_shift)
theta_third_truth <- function(a, p = dgp_parameters()) -cos(a + p$dose_shift)
theta_fourth_truth <- function(a, p = dgp_parameters()) sin(a + p$dose_shift)
theta_fifth_truth <- function(a, p = dgp_parameters()) cos(a + p$dose_shift)

exposure_marginal_variance <- function(p = dgp_parameters()) {
  p$sigma_x^2 * sum(p$alpha_x^2) + p$alpha_u^2 + p$sigma_a^2
}

orthogonal_dose_constant <- function(p = dgp_parameters()) {
  exp(-exposure_marginal_variance(p) / 2)
}

outcome_dose_basis <- function(a, p = dgp_parameters()) {
  constant <- orthogonal_dose_constant(p)
  cbind(A = a, S_A = sin(a) - constant * a, C_A = cos(a) - constant)
}

outcome_bridge_truth <- function(w, a, x1, x2, p = dgp_parameters()) {
  theta_truth(a, p) + tanh(w)
}

simulate_proximal_data <- function(n, seed, keep_latent = FALSE) {
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
    z <- p$zeta0 + p$zeta_a * a + p$zeta_x[[1L]] * x1 +
      p$zeta_x[[2L]] * x2 + p$zeta_u * u + eps_z
    w_mean <- p$mu0 + p$mu_x[[1L]] * x1 + p$mu_x[[2L]] * x2 + p$mu_u * u
    w <- w_mean + eps_w
    y <- theta_truth(a, p) + tanh(w) + eps_y
    data.frame(Y = y, A = a, Z = z, W = w, X1 = x1, X2 = x2,
               U = u, W_mean = w_mean, eps_Y = eps_y,
               stringsAsFactors = FALSE)
  })
  if (!keep_latent) out[c("U", "W_mean", "eps_Y")] <- NULL
  out
}

treatment_bridge_moment_threshold <- function(p = dgp_parameters()) {
  variance_a <- exposure_marginal_variance(p)
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
  variance_a <- exposure_marginal_variance(p)
  denominator <- p$sigma_a^2 * p$zeta_u^2 + p$alpha_u^2 * p$sigma_z^2
  lambda <- p$alpha_u^2 / (2 * denominator)
  gaussian_integrability <- p$sigma_a^2 * p$zeta_u^2 / denominator
  mean_a_x <- p$alpha0 + p$alpha_x[[1L]] * x1 + p$alpha_x[[2L]] * x2
  residual <- z - p$zeta0 - p$zeta_a * a - p$zeta_x[[1L]] * x1 -
    p$zeta_x[[2L]] * x2 - (p$zeta_u / p$alpha_u) * (a - mean_a_x)
  log_q <- log(p$sigma_a / sqrt(variance_a)) +
    0.5 * log(gaussian_integrability) - a^2 / (2 * variance_a) +
    lambda * residual^2
  exp(log_q)
}

dgp_diagnostics <- function(data, target_range = c(-0.9, 0.9)) {
  p <- dgp_parameters()
  matrix_data <- as.matrix(data[c("Y", "A", "Z", "W", "X1", "X2")])
  c(
    n = nrow(data),
    marginal_variance_a = exposure_marginal_variance(p),
    q_moment_threshold = treatment_bridge_moment_threshold(p),
    mean_a = mean(data$A), sd_a = sd(data$A),
    q01_a = unname(quantile(data$A, 0.01)),
    q99_a = unname(quantile(data$A, 0.99)),
    fraction_in_target_range = mean(data$A >= target_range[[1L]] &
                                      data$A <= target_range[[2L]]),
    max_abs_y = max(abs(data$Y)),
    y_bound_satisfied = as.numeric(max(abs(data$Y)) <= 5 + 1e-12),
    min_abs_theta3_grid = min(abs(theta_third_truth(seq(target_range[[1L]],
                                                        target_range[[2L]], length.out = 1001L)))),
    all_finite = as.numeric(all(is.finite(matrix_data)))
  )
}
