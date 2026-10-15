root <- Sys.getenv("PDR_JOINT_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- commandArgs(FALSE)[startsWith(commandArgs(FALSE), "--file=")][[1L]]
  root <- dirname(dirname(sub("^--file=", "", file_arg)))
}
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "R", "utils.R"), local = .GlobalEnv)
source_simulation_files(root, compile_cpp = TRUE)
cli <- parse_cli()
quick <- cli_value(cli, "quick", FALSE, "logical")
cfg <- simulation_config("test")

near <- function(x, y, tolerance, label) {
  error <- max(abs(x - y))
  if (!is.finite(error) || error > tolerance)
    stop(sprintf("%s failed: max error %.3e > %.3e", label, error, tolerance), call. = FALSE)
  message(sprintf("PASS %-52s max error %.3e", label, error))
}

near_relative <- function(x, y, tolerance, label) {
  error <- max(abs(x - y) / pmax(1, abs(y)))
  if (!is.finite(error) || error > tolerance)
    stop(sprintf("%s failed: max scaled error %.3e > %.3e", label, error, tolerance),
         call. = FALSE)
  message(sprintf("PASS %-52s max scaled error %.3e", label, error))
}

# Analytic truth, bounded outcome, and active derivative bias.
grid <- seq(-0.9, 0.9, length.out = 1001L)
step <- 1e-5
numeric_derivative <- (theta_truth(grid + step) - theta_truth(grid - step)) / (2 * step)
near(numeric_derivative, theta_derivative_truth(grid), 2e-10,
     "analytic first derivative")
assert_true(min(abs(theta_third_truth(grid))) > 0.45,
            "theta''' is not active throughout the reporting grid")
latent <- simulate_proximal_data(if (quick) 20000L else 100000L, 891L, keep_latent = TRUE)
assert_true(max(abs(latent$Y)) <= 5 + 1e-12, "bounded-outcome construction failed")
near(latent$Y - outcome_bridge_truth(latent$W, latent$A, latent$X1, latent$X2),
     latent$eps_Y, 1e-13, "outcome bridge residual identity")

# Exact representation of theta in the stabilized outcome basis.
p <- dgp_parameters(); constant <- orthogonal_dose_constant(p)
coefficient <- c(`(Intercept)` = 2.5 + constant * sin(p$dose_shift),
                 A = constant * cos(p$dose_shift), X1 = 0, X2 = 0,
                 S_A = cos(p$dose_shift), C_A = sin(p$dose_shift))
base <- outcome_base_design(data.frame(A = grid, X1 = 0, X2 = 0, W = 0, Z = 0))
near(drop(base %*% coefficient), theta_truth(grid), 2e-15,
     "outcome dose-basis representation")

# The exact treatment bridge remains in the complete raw quadratic log model.
sample <- latent[seq_len(min(5000L, nrow(latent))), ]
q <- treatment_bridge_truth(sample$Z, sample$A, sample$X1, sample$X2)
basis <- quadratic_basis(sample$Z, sample$A, sample$X1, sample$X2)
fit <- lm.fit(basis, log(q))
near(drop(basis %*% fit$coefficients), log(q), 2e-10,
     "log-quadratic treatment bridge representation")
assert_true(treatment_bridge_moment_threshold() > 4,
            "Treatment bridge lacks the fourth moment required by diagnostics")

# Strong latent treatment-bridge equation at a fixed (u,a,x).
u0 <- 0.4; a0 <- 0.3; x10 <- -0.1; x20 <- 0.15
mc_n <- if (quick) 100000L else 500000L
z_mc <- with_preserved_seed(733L, {
  rnorm(mc_n, p$zeta0 + p$zeta_a * a0 + p$zeta_x[[1L]] * x10 +
          p$zeta_x[[2L]] * x20 + p$zeta_u * u0, p$sigma_z)
})
lhs <- mean(treatment_bridge_truth(z_mc, a0, x10, x20))
mean_ax <- p$alpha0 + p$alpha_x[[1L]] * x10 + p$alpha_x[[2L]] * x20
rhs <- dnorm(a0, 0, sqrt(exposure_marginal_variance(p))) /
  dnorm(a0, mean_ax + p$alpha_u * u0, p$sigma_a)
near(lhs, rhs, if (quick) 8e-3 else 4e-3, "latent treatment-bridge equation")

# Population outcome moment system is nonsingular and omit-W is genuinely misspecified.
instrument <- outcome_instruments(sample)
correct_design <- outcome_design(sample, "correct")
singular_values <- svd(crossprod(instrument, correct_design) / nrow(sample),
                       nu = 0, nv = 0)$d
assert_true(min(singular_values) > 1e-3, "Correct outcome moment system is nearly singular")
omit_design <- outcome_design(sample, "omit_w")
omega <- crossprod(instrument) / nrow(sample)
g <- crossprod(instrument, omit_design) / nrow(sample)
gy <- crossprod(instrument, sample$Y) / nrow(sample)
beta_m <- solve(crossprod(g, solve(omega, g)), crossprod(g, solve(omega, gy)))
omit_moment <- drop(crossprod(instrument, sample$Y - drop(omit_design %*% beta_m)) /
                       nrow(sample))
assert_true(max(abs(omit_moment)) > 1e-3, "Omit-W outcome class is not empirically distinct")

# Direct b=h identities for the level and derivative robust-bias corrections.
set.seed(11)
a <- sort(rnorm(600)); y <- 1 + 0.4 * a - 0.3 * a^2 + 0.2 * a^3 + sin(a)
eval <- 0.17; h <- 0.64; n <- length(a); u <- (a - eval) / h
kh <- kernel_epanechnikov(u) / h
w1 <- cbind(1, u); w2 <- cbind(1, u, u^2); w3 <- cbind(1, u, u^2, u^3)
d1 <- crossprod(w1, w1 * kh) / n; d2 <- crossprod(w2, w2 * kh) / n
d3 <- crossprod(w3, w3 * kh) / n
i1 <- solve(d1); i2 <- solve(d2); i3 <- solve(d3)
r0 <- drop(w1 %*% i1[, 1L]) * kh
c2 <- drop(i1[1L, , drop = FALSE] %*% colMeans(w1 * (kh * u^2)))
t2 <- drop(w2 %*% i2[, 3L]) * kh
explicit_level <- r0 - c2 * t2
direct_level <- drop(w2 %*% i2[, 1L]) * kh
near(explicit_level, direct_level, 2e-12, "b=h level RBC equals local-quadratic intercept")

r1 <- drop(w2 %*% i2[, 2L]) * kh / h
v3 <- drop(i2 %*% colMeans(w2 * (kh * u^3)))
c3 <- v3[[2L]]
s3 <- 6 * drop(w3 %*% i3[, 4L]) * kh / h^3
explicit_derivative <- r1 - h^2 * c3 * s3 / 6
direct_derivative <- drop(w3 %*% i3[, 2L]) * kh / h
near(explicit_derivative, direct_derivative, 3e-12,
     "b=h derivative RBC equals local-cubic slope")
m2 <- colMeans(w2 * (kh * y)); m3 <- colMeans(w3 * (kh * y))
beta2 <- drop(i2 %*% m2); beta3 <- drop(i3 %*% m3)
gamma1 <- r1 * drop(w2 %*% beta2)
gamma3 <- s3 * drop(w3 %*% beta3)
gamma_c <- drop(w2 %*% i2[, 2L]) * kh * (u^3 - drop(w2 %*% v3))
theta3 <- 6 * beta3[[4L]] / h^3
explicit_gamma <- gamma1 - h^2 * c3 * gamma3 / 6 + h^2 * theta3 * gamma_c / 6
direct_gamma <- direct_derivative * drop(w3 %*% beta3)
near(explicit_gamma, direct_gamma, 2e-11,
     "three-term derivative design correction identity")
component <- local_polynomial_components(a, y, eval, h, 3L, 1L, cfg)
near(component$Gamma, direct_derivative, 2e-12, "production local-cubic derivative weight")
near(component$reproduction, c(0, 1, 0, 0), 2e-12,
     "cubic derivative polynomial reproduction")

# Same-fold product-measure term as an empirical Gateaux derivative.  For an
# additive bridge h(w,a,x)=g(a)+r(w,x), perturbing the empirical law of (W,X)
# within a fold produces mean(Gamma) * {r_i-mean(r)}.
set.seed(312)
fold_check <- rep(seq_len(cfg$folds), length.out = n)
r_check <- tanh(rnorm(n)) + 0.2 * rnorm(n)
analytic_product <- numeric(n); numeric_product <- numeric(n); epsilon <- 1e-7
for (fold_id in seq_len(cfg$folds)) {
  index <- which(fold_check == fold_id)
  mean_gamma <- mean(component$Gamma[index]); mean_r <- mean(r_check[index])
  analytic_product[index] <- mean_gamma * (r_check[index] - mean_r)
  numeric_product[index] <- vapply(index, function(i) {
    perturbed <- (1 - epsilon) * mean_r + epsilon * r_check[[i]]
    (mean_gamma * perturbed - mean_gamma * mean_r) / epsilon
  }, numeric(1L))
}
near(analytic_product, numeric_product, 2e-9,
     "same-fold product-measure Gateaux derivative")

# Epanechnikov equivalent-kernel constants for the primary and RBC slopes.
integration <- seq(-1, 1, length.out = 200001L)
k <- kernel_epanechnikov(integration)
l2 <- 5 * integration * k
l3 <- (175 * integration - 315 * integration^3) * k / 8
near(trapz_grid(integration, l2^2), 15 / 7, 2e-9,
     "degree-two slope equivalent-kernel L2 norm")
near(trapz_grid(integration, l3^2), 525 / 44, 2e-8,
     "degree-three slope equivalent-kernel L2 norm")
near(vapply(0:3, function(power) trapz_grid(integration, l3 * integration^power),
             numeric(1L)), c(0, 1, 0, 0), 1e-9,
     "degree-three equivalent-kernel moments")

# Accelerated level LOOCV against literal leave-one-out refits.
set.seed(94); a_cv <- sort(rnorm(140)); xi_cv <- cbind(CC = sin(a_cv) + rnorm(140, sd=.3),
                                                       MC = cos(a_cv) + rnorm(140, sd=.4))
cfg_cv <- cfg; cfg_cv$cv_min_neighbors <- 5L; cfg_cv$level_cv_trim <- 0.1
bandwidths <- c(0.45, 0.7)
fast_level <- level_loocv(a_cv, xi_cv, bandwidths, cfg_cv)
trim <- quantile(a_cv, c(.1,.9), names=FALSE)
brute_level <- function(h, degree, column) {
  residual <- numeric()
  for (i in which(a_cv >= trim[[1L]] & a_cv <= trim[[2L]])) {
    keep <- setdiff(seq_along(a_cv), i); uu <- (a_cv[keep] - a_cv[[i]]) / h
    weight <- kernel_epanechnikov(uu); active <- weight > 0
    if (sum(active) < cfg_cv$cv_min_neighbors) next
    ww <- outer(uu[active], 0:degree, `^`)
    beta <- solve(crossprod(ww, ww * weight[active]),
                  crossprod(ww, weight[active] * xi_cv[keep[active], column]))
    residual <- c(residual, xi_cv[i, column] - beta[[1L]])
  }
  mean(residual^2)
}
literal_ll <- sapply(seq_len(ncol(xi_cv)), function(j)
  vapply(bandwidths, function(value) brute_level(value, 1L, j), numeric(1L)))
literal_lq <- sapply(seq_len(ncol(xi_cv)), function(j)
  vapply(bandwidths, function(value) brute_level(value, 2L, j), numeric(1L)))
near(fast_level$risk_ll, literal_ll, 5e-12, "accelerated local-linear level LOOCV")
near(fast_level$risk_pdr, literal_lq, 5e-12, "accelerated quadratic-intercept level LOOCV")

# Accelerated derivative CV against the same analytic center derivative in R.
density <- derivative_density(a_cv, cfg_cv)
fast_derivative <- derivative_loocv(a_cv, xi_cv, density$inverse, bandwidths, cfg_cv)
taper_r <- function(a) {
  if (a <= -1 || a >= 1) return(c(rho=0,rhop=0))
  if (a >= -.9 && a <= .9) return(c(rho=1,rhop=0))
  if (a < -.9) { t <- (-.9-a)/.1; direction <- -1 } else { t <- (a-.9)/.1; direction <- 1 }
  c(rho=1-10*t^3+15*t^4-6*t^5,
    rhop=(-30*t^2+60*t^3-30*t^4)*direction/.1)
}
brute_derivative <- function(h, degree, column) {
  value <- numeric()
  for (i in which(a_cv > -1 & a_cv < 1)) {
    keep <- setdiff(seq_along(a_cv), i); uu <- (a_cv[keep] - a_cv[[i]]) / h
    weight <- kernel_epanechnikov(uu); active <- weight > 0
    if (sum(active) < cfg_cv$cv_min_neighbors) next
    uu <- uu[active]; weight <- weight[active]; yy <- xi_cv[keep[active], column]
    up <- vapply(0:(2*degree+2L), function(power) sum(uu^power), numeric(1L))
    t <- vapply(0:(2*degree), function(j) up[[j+1L]]-up[[j+3L]], numeric(1L))
    td <- vapply(0:(2*degree), function(j) {
      first <- if (j == 0) 0 else j * up[[j]]
      -(first - (j+2) * up[[j+2L]]) / h
    }, numeric(1L))
    yup <- vapply(0:(degree+2L), function(power) sum(yy*uu^power), numeric(1L))
    m <- vapply(0:degree, function(j) yup[[j+1L]] - yup[[j+3L]], numeric(1L))
    d <- outer(0:degree, 0:degree, function(r,c) t[r+c+1L])
    beta <- drop(solve(d, m))
    md <- vapply(0:degree, function(j) {
      first <- if (j == 0) 0 else j * yup[[j]]
      -(first - (j+2) * yup[[j+2L]]) / h
    }, numeric(1L))
    dd <- outer(0:degree,0:degree,function(r,c) td[r+c+1L])
    beta_dot <- solve(d, md - dd %*% beta)
    g <- beta[[2L]]/h; gd <- beta_dot[[2L]]/h; taper <- taper_r(a_cv[[i]])
    value <- c(value, density$inverse[[i]] *
                 (taper[["rho"]]*g^2 + 2*xi_cv[i,column]*
                    (taper[["rho"]]*gd+taper[["rhop"]]*g)))
  }
  mean(value)
}
literal_dll <- sapply(seq_len(ncol(xi_cv)), function(j)
  vapply(bandwidths, function(value) brute_derivative(value, 1L, j), numeric(1L)))
literal_dpdr <- sapply(seq_len(ncol(xi_cv)), function(j)
  vapply(bandwidths, function(value) brute_derivative(value, 3L, j), numeric(1L)))
near(fast_derivative$risk_ll, literal_dll, 2e-10,
     "accelerated derivative-CV local-linear slope")
near_relative(fast_derivative$risk_pdr, literal_dpdr, 5e-12,
              "accelerated derivative-CV debiased local-cubic slope")

# Full seed-reproducible pipeline, covariance scale, and artifact validator.
cfg$gaussian_band_draws <- 120L
one <- run_one_replication(800L, 1L, cfg, root)
two <- run_one_replication(800L, 1L, cfg, root)
assert_true(identical(one$curve, two$curve) && identical(one$band, two$band) &&
              identical(one$bandwidth, two$bandwidth),
            "Repeated seed-matched runs are not identical")
assert_true(max(one$curve$reproduction_error) < 1e-10,
            "Full pipeline violates local-polynomial reproduction")
temporary <- tempfile(fileext = ".rds"); on.exit(unlink(temporary), add = TRUE)
saveRDS(one, temporary)
assert_true(validate_part(temporary, cfg, 800L, 1L), "Full artifact validation failed")
message("PASS full joint pipeline, deterministic seeds, covariance, and schema")
message("All mathematical and implementation validation checks passed.")
