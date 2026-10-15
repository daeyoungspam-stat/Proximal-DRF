parent_root <- normalizePath(Sys.getenv("PDR_JOINT_ROOT", unset = ""), mustWork = TRUE)
source(file.path(parent_root, "R", "utils.R"), local = .GlobalEnv)
source(file.path(parent_root, "proxlq_addon", "R", "addon.R"), local = .GlobalEnv)
source_proxlq_addon(parent_root, compile_cpp = TRUE)
cli <- parse_cli()
quick <- cli_value(cli, "quick", FALSE, "logical")

near <- function(x, y, tolerance, label) {
  difference <- max(abs(x - y))
  if (!is.finite(difference) || difference > tolerance)
    stop(sprintf("%s failed: max absolute difference %.3e > %.3e",
                 label, difference, tolerance), call. = FALSE)
  message(sprintf("PASS %s (%.3e)", label, difference))
}

near_relative <- function(x, y, tolerance, label) {
  scale <- pmax(1, abs(x), abs(y))
  difference <- max(abs(x - y) / scale)
  if (!is.finite(difference) || difference > tolerance)
    stop(sprintf("%s failed: max relative difference %.3e > %.3e",
                 label, difference, tolerance), call. = FALSE)
  message(sprintf("PASS %s (%.3e)", label, difference))
}

cfg <- proxlq_addon_config("test", parent_root)
validate_proxlq_design(cfg)
assert_true(identical(cfg$parent_config_id, config_fingerprint(simulation_config("production"), parent_root)),
            "Frozen parent fingerprint changed")
message("PASS frozen parent fingerprint")

# Degree-two slope reproduction and exact differentiation of a quadratic.
with_preserved_seed(2026082301L, {
  a <- rnorm(1200)
  eval <- 0.13
  y <- 1.7 + 2.2 * a - 0.8 * a^2
  component <- local_polynomial_components(
    a, y, eval, 0.55, degree = 2L, derivative_order = 1L, cfg$base_cfg
  )
  near(component$reproduction, c(0, 1, 0), 2e-12,
       "degree-two derivative polynomial reproduction")
  near(component$theta, 2.2 - 1.6 * eval, 2e-12,
       "noise-free quadratic derivative")
})

# Epanechnikov degree-two derivative equivalent-kernel constants.
epan <- function(u) 0.75 * (1 - u^2) * (abs(u) <= 1)
l2 <- function(u) 5 * u * epan(u)
near(integrate(function(u) l2(u), -1, 1)$value, 0, 2e-12,
     "equivalent-kernel constant moment")
near(integrate(function(u) l2(u) * u, -1, 1)$value, 1, 2e-12,
     "equivalent-kernel linear moment")
near(integrate(function(u) l2(u) * u^2, -1, 1)$value, 0, 2e-12,
     "equivalent-kernel quadratic moment")
near(integrate(function(u) l2(u)^2, -1, 1)$value, 15 / 7, 2e-10,
     "equivalent-kernel L2 norm")
near(integrate(function(u) l2(u) * u^3, -1, 1)$value, 3 / 7, 2e-10,
     "equivalent-kernel cubic moment")

# Independent integration-by-parts verification of the derivative-risk sign.
taper <- function(a) {
  rho <- rhop <- numeric(length(a))
  flat <- a >= -0.9 & a <= 0.9
  rho[flat] <- 1
  left <- a > -1 & a < -0.9
  right <- a > 0.9 & a < 1
  for (index in which(left | right)) {
    if (left[[index]]) {
      t <- (-0.9 - a[[index]]) / 0.1; direction <- -1
    } else {
      t <- (a[[index]] - 0.9) / 0.1; direction <- 1
    }
    rho[[index]] <- 1 - 10 * t^3 + 15 * t^4 - 6 * t^5
    rhop[[index]] <- (-30 * t^2 + 60 * t^3 - 30 * t^4) * direction / 0.1
  }
  list(rho = rho, rhop = rhop)
}
theta <- function(a) 2.5 + sin(a + 0.2)
theta_prime <- function(a) cos(a + 0.2)
g <- function(a) 0.4 + 0.7 * a - 0.3 * a^2
g_prime <- function(a) 0.7 - 0.6 * a
integrate_pieces <- function(fun) sum(vapply(
  list(c(-1, -0.9), c(-0.9, 0.9), c(0.9, 1)),
  function(bounds) integrate(fun, bounds[[1L]], bounds[[2L]],
                             subdivisions = 1000L, rel.tol = 1e-12,
                             abs.tol = 1e-12)$value,
  numeric(1L)
))
q_left <- integrate_pieces(function(a) {
  tp <- taper(a)
  tp$rho * g(a)^2 + 2 * theta(a) * (tp$rhop * g(a) + tp$rho * g_prime(a))
})
q_right <- integrate_pieces(function(a) {
  tp <- taper(a)
  tp$rho * ((g(a) - theta_prime(a))^2 - theta_prime(a)^2)
})
near(q_left, q_right, 2e-11, "integration-by-parts derivative-risk identity")
boundary <- taper(c(-1, 1))
near(c(boundary$rho, boundary$rhop), rep(0, 4), 1e-14,
     "taper boundary conditions")

# Accelerated p=2 criterion against literal leave-one-out refits.
cfg_cv <- proxlq_addon_config("test", parent_root)
cfg_cv$base_cfg$cv_min_neighbors <- 8L
with_preserved_seed(2026082302L, {
  a_cv <- rnorm(260, sd = 0.72)
  xi_cv <- cbind(
    first = 2.5 + sin(a_cv + 0.2) + rnorm(length(a_cv), sd = 0.4),
    second = 1.2 - 0.4 * a_cv + 0.6 * a_cv^2 + rnorm(length(a_cv), sd = 0.6)
  )
  bandwidths <- c(0.31, 0.39, 0.49)
  density <- derivative_density(a_cv, cfg_cv$base_cfg)
  accelerated <- derivative_lq_loocv(
    a_cv, xi_cv, density$inverse, bandwidths, cfg_cv
  )

  brute <- function(h, column) {
    value <- numeric()
    for (i in which(a_cv > -1 & a_cv < 1)) {
      keep <- setdiff(seq_along(a_cv), i)
      u <- (a_cv[keep] - a_cv[[i]]) / h
      weight <- kernel_epanechnikov(u)
      active <- weight > 0
      if (sum(active) < cfg_cv$base_cfg$cv_min_neighbors) next
      u <- u[active]; weight <- weight[active]
      response <- xi_cv[keep[active], column]
      up <- vapply(0:6, function(power) sum(u^power), numeric(1L))
      t <- vapply(0:4, function(j) up[[j + 1L]] - up[[j + 3L]], numeric(1L))
      tdot <- vapply(0:4, function(j) {
        first <- if (j == 0) 0 else j * up[[j]]
        -(first - (j + 2) * up[[j + 2L]]) / h
      }, numeric(1L))
      yup <- vapply(0:4, function(power) sum(response * u^power), numeric(1L))
      moment <- vapply(0:2, function(j) yup[[j + 1L]] - yup[[j + 3L]], numeric(1L))
      d <- outer(0:2, 0:2, function(r, c) t[r + c + 1L])
      beta <- drop(solve(d, moment))
      moment_dot <- vapply(0:2, function(j) {
        first <- if (j == 0) 0 else j * yup[[j]]
        -(first - (j + 2) * yup[[j + 2L]]) / h
      }, numeric(1L))
      d_dot <- outer(0:2, 0:2, function(r, c) tdot[r + c + 1L])
      beta_dot <- drop(solve(d, moment_dot - d_dot %*% beta))
      derivative <- beta[[2L]] / h
      derivative_dot <- beta_dot[[2L]] / h
      tp <- taper(a_cv[[i]])
      value <- c(value, density$inverse[[i]] *
        (tp$rho * derivative^2 + 2 * xi_cv[i, column] *
         (tp$rho * derivative_dot + tp$rhop * derivative)))
    }
    mean(value)
  }
  literal <- sapply(seq_len(ncol(xi_cv)), function(column)
    vapply(bandwidths, brute, numeric(1L), column = column))
  near_relative(accelerated$risk_lq, literal, 3e-11,
                "accelerated degree-two derivative risk")

  permutation <- sample.int(length(a_cv))
  permuted <- derivative_lq_loocv(
    a_cv[permutation], xi_cv[permutation, , drop = FALSE],
    density$inverse[permutation], bandwidths, cfg_cv
  )
  near_relative(accelerated$risk_lq, permuted$risk_lq, 2e-12,
                "degree-two CV permutation invariance")
})

# Exact undersmoothing formulas and common raw selector.
for (n in c(2500L, 5000L, 7500L, 10000L)) {
  raw <- 1.37 * n^(-1 / 7)
  h10 <- raw / proxlq_divisor(n, "ProxLQ-US-log10")
  h27 <- raw / proxlq_divisor(n, "ProxLQ-US-log2over7")
  near(n * h10^7, 1.37^7 / log10(n)^7, 2e-14,
       paste0("log10 nh7 identity n=", n))
  near(n * h27^7, 1.37^7 / log(n)^2, 2e-14,
       paste0("log-power nh7 identity n=", n))
  assert_true(0 < h10 && h10 < h27 && h27 < raw,
              "Undersmoothing bandwidth ordering failed")
}
message("PASS exact undersmoothing transformations")

# Selector mechanics: finite-candidate endpoint flags and deterministic tie break.
selection <- select_from_allowed(
  c(0.1, 0.2, 0.3), matrix(c(3, 1, 1), ncol = 1L,
                            dimnames = list(NULL, "CC")), list(c(0.1, 0.2, 0.3))
)
near(selection$h[[1L]], 0.2, 1e-14, "CV tie break to smaller bandwidth")
assert_true(!selection$lower_hit[[1L]] && !selection$upper_hit[[1L]],
            "Interior selection boundary flags failed")

# End-to-end deterministic schema test.
test_cfg <- proxlq_addon_config("test", parent_root)
one <- run_one_proxlq_replication(800L, 1L, test_cfg, parent_root)
two <- run_one_proxlq_replication(800L, 1L, test_cfg, parent_root)
assert_true(identical(one$curve, two$curve) &&
              identical(one$bandwidth, two$bandwidth) &&
              identical(one$cv_diagnostics, two$cv_diagnostics),
            "Same-seed ProxLQ runs are not identical")
temporary <- tempfile(fileext = ".rds")
on.exit(unlink(temporary), add = TRUE)
saveRDS(one, temporary)
assert_true(validate_proxlq_part(temporary, test_cfg, 800L, 1L, parent_root),
            "End-to-end ProxLQ artifact validation failed")
message("PASS deterministic end-to-end ProxLQ pipeline")

if (!quick) {
  production_cfg <- proxlq_addon_config("production", parent_root)
  matched <- run_one_proxlq_replication(2500L, 1L, production_cfg, parent_root)
  base_path <- file.path(parent_root, "results", "config_4e7d28a1", "n_2500",
                         "rep_0001.rds")
  assert_true(file.exists(base_path), "Frozen base artifact is missing")
  base <- readRDS(base_path)
  assert_true(identical(matched$metadata$dgp_seed, base$metadata$dgp_seed) &&
                identical(matched$metadata$fold_seed, base$metadata$fold_seed),
              "Matched production seeds differ from frozen base")
  near(unlist(matched$dgp_diagnostics), unlist(base$dgp_diagnostics), 0,
       "matched DGP diagnostics")
  numeric_columns <- names(matched$nuisance_diagnostics)[
    vapply(matched$nuisance_diagnostics, is.numeric, logical(1L))]
  nuisance_new <- as.matrix(matched$nuisance_diagnostics[numeric_columns])
  nuisance_old <- as.matrix(base$nuisance_diagnostics[numeric_columns])
  assert_true(identical(is.na(nuisance_new), is.na(nuisance_old)),
              "Matched nuisance NA patterns differ")
  near(nuisance_new[!is.na(nuisance_new)], nuisance_old[!is.na(nuisance_old)], 0,
       "matched nuisance diagnostics")
  message("PASS frozen-base matched replication")
}

message("All ProxLQ add-on mathematical and implementation validations passed.")
