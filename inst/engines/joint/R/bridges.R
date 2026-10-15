quadratic_basis <- function(proxy, a, x1, x2) {
  cbind(
    `(Intercept)` = 1,
    proxy = proxy, A = a, X1 = x1, X2 = x2,
    proxy2 = proxy^2, A2 = a^2, X1_2 = x1^2, X2_2 = x2^2,
    proxy_A = proxy * a,
    proxy_X1 = proxy * x1, proxy_X2 = proxy * x2,
    A_X1 = a * x1, A_X2 = a * x2, X1_X2 = x1 * x2
  )
}

main_effects_basis <- function(proxy, a, x1, x2) {
  cbind(`(Intercept)` = 1, proxy = proxy, A = a, X1 = x1, X2 = x2)
}

treatment_basis <- function(data, model = c("quadratic", "main_effects"),
                            proxy = c("Z", "W")) {
  model <- match.arg(model)
  proxy <- match.arg(proxy)
  if (model == "quadratic")
    quadratic_basis(data[[proxy]], data$A, data$X1, data$X2)
  else
    main_effects_basis(data[[proxy]], data$A, data$X1, data$X2)
}

product_calibration_target <- function(data) {
  mw <- mean(data$W); ma <- mean(data$A)
  mx1 <- mean(data$X1); mx2 <- mean(data$X2)
  c(
    `(Intercept)` = 1,
    proxy = mw, A = ma, X1 = mx1, X2 = mx2,
    proxy2 = mean(data$W^2), A2 = mean(data$A^2),
    X1_2 = mean(data$X1^2), X2_2 = mean(data$X2^2),
    proxy_A = mw * ma,
    proxy_X1 = mean(data$W * data$X1), proxy_X2 = mean(data$W * data$X2),
    A_X1 = ma * mx1, A_X2 = ma * mx2, X1_X2 = mean(data$X1 * data$X2)
  )
}

outcome_base_design <- function(data, a = data$A,
                                x1 = data$X1, x2 = data$X2) {
  dose <- outcome_dose_basis(a)
  cbind(`(Intercept)` = 1, A = dose[, "A"], X1 = x1, X2 = x2,
        S_A = dose[, "S_A"], C_A = dose[, "C_A"])
}

outcome_design <- function(data, model = c("correct", "omit_w"), a = data$A,
                           w = data$W, x1 = data$X1, x2 = data$X2) {
  model <- match.arg(model)
  base <- outcome_base_design(data, a, x1, x2)
  if (model == "correct") cbind(base, tanh_W = tanh(w)) else base
}

outcome_instruments <- function(data) {
  cbind(outcome_base_design(data), Z = data$Z)
}

fit_outcome_bridge <- function(data, model = c("correct", "omit_w"), cfg) {
  model <- match.arg(model)
  hmat <- outcome_design(data, model)
  zmat <- outcome_instruments(data)
  omega <- crossprod(zmat) / nrow(data)
  omega_inverse <- safe_solve(
    omega, ridge = cfg$nuisance_ridge,
    condition_limit = cfg$nuisance_condition_limit
  )
  g <- crossprod(zmat, hmat) / nrow(data)
  gy <- crossprod(zmat, data$Y) / nrow(data)
  lhs <- crossprod(g, omega_inverse %*% g)
  rhs <- crossprod(g, omega_inverse %*% gy)
  coefficient <- drop(safe_solve(
    lhs, rhs, ridge = cfg$nuisance_ridge,
    condition_limit = cfg$nuisance_condition_limit
  ))
  names(coefficient) <- colnames(hmat)
  moment <- drop(crossprod(zmat, data$Y - drop(hmat %*% coefficient)) / nrow(data))
  score <- drop(crossprod(g, omega_inverse %*% moment))
  list(
    model = model,
    coefficients = coefficient,
    max_abs_moment = max(abs(moment)),
    max_abs_score = max(abs(score)),
    moments = setNames(moment, colnames(zmat)),
    condition = kappa(lhs, exact = TRUE)
  )
}

predict_outcome_bridge <- function(fit, data, a = data$A) {
  drop(outcome_design(data, fit$model, a = a) %*% fit$coefficients)
}

standardize_columns <- function(x) {
  scale <- sqrt(colMeans(x^2))
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  scale[[1L]] <- 1
  list(x = sweep(x, 2L, scale, "/"), scale = scale)
}

fit_treatment_bridge <- function(data, model = c("quadratic", "main_effects"), cfg) {
  model <- match.arg(model)
  b_raw <- treatment_basis(data, model, proxy = "Z")
  g_raw <- treatment_basis(data, "quadratic", proxy = "W")
  target_raw <- product_calibration_target(data)
  assert_true(identical(colnames(g_raw), names(target_raw)),
              "Treatment calibration target/basis mismatch")

  b_standard <- standardize_columns(b_raw)
  g_standard <- standardize_columns(g_raw)
  b <- b_standard$x
  g <- g_standard$x
  target <- target_raw / g_standard$scale
  n <- nrow(data); parameter_count <- ncol(b); moment_count <- ncol(g)

  evaluate <- function(parameter, jacobian = TRUE) {
    eta <- drop(b %*% parameter)
    if (any(!is.finite(eta)) || max(abs(eta)) > cfg$q_eta_limit) {
      return(list(moment = rep(Inf, moment_count),
                  jacobian = matrix(NA_real_, moment_count, parameter_count),
                  eta = eta, q = rep(NA_real_, n)))
    }
    q <- exp(eta)
    moment <- colMeans(g * q) - target
    jac <- if (jacobian) crossprod(g * q, b) / n else NULL
    list(moment = moment, jacobian = jac, eta = eta, q = q)
  }

  parameter <- rep(0, parameter_count)
  iterations <- 0L; line_search_failures <- 0L
  for (iteration in seq_len(cfg$q_max_iter)) {
    iterations <- iteration
    current <- evaluate(parameter)
    norm0 <- sqrt(sum(current$moment^2))
    score0 <- if (all(is.finite(current$jacobian)))
      drop(crossprod(current$jacobian, current$moment)) else rep(Inf, parameter_count)
    exact <- moment_count == parameter_count
    if (is.finite(norm0) &&
        ((exact && max(abs(current$moment)) <= cfg$q_root_tolerance) ||
         (!exact && max(abs(score0)) <= cfg$q_root_tolerance))) break
    if (!all(is.finite(current$jacobian))) break
    step <- tryCatch(
      if (exact) {
        drop(qr.solve(current$jacobian, current$moment, tol = 1e-12))
      } else {
        normal <- crossprod(current$jacobian)
        drop(qr.solve(normal + diag(cfg$nuisance_ridge, parameter_count),
                      crossprod(current$jacobian, current$moment), tol = 1e-12))
      },
      error = function(e) rep(NA_real_, parameter_count)
    )
    if (!all(is.finite(step))) break
    accepted <- FALSE; step_scale <- 1
    for (line in 0:24) {
      proposal <- parameter - step_scale * step
      trial <- evaluate(proposal, jacobian = FALSE)
      norm1 <- sqrt(sum(trial$moment^2))
      if (is.finite(norm1) && norm1 < norm0 * (1 - 1e-4 * step_scale)) {
        parameter <- proposal; accepted <- TRUE; break
      }
      step_scale <- step_scale / 2
    }
    if (!accepted) { line_search_failures <- line_search_failures + 1L; break }
  }

  final <- evaluate(parameter)
  final_score <- if (all(is.finite(final$jacobian)))
    drop(crossprod(final$jacobian, final$moment)) else rep(Inf, parameter_count)
  exact <- moment_count == parameter_count
  needs_fallback <- !all(is.finite(final$moment)) ||
    (exact && max(abs(final$moment)) > cfg$q_accept_tolerance) ||
    (!exact && max(abs(final_score)) > cfg$q_accept_tolerance)
  if (needs_fallback) {
    objective <- function(value) {
      result <- evaluate(value)
      if (any(!is.finite(result$moment))) return(1e100)
      sum(result$moment^2)
    }
    gradient <- function(value) {
      result <- evaluate(value)
      if (any(!is.finite(result$moment)) || any(!is.finite(result$jacobian)))
        return(rep(0, parameter_count))
      drop(2 * crossprod(result$jacobian, result$moment))
    }
    optimization <- optim(parameter, objective, gradient, method = "BFGS",
                          control = list(maxit = 500L, reltol = 1e-12))
    if (is.finite(optimization$value) && optimization$value < objective(parameter))
      parameter <- optimization$par
    final <- evaluate(parameter)
  }

  final_score <- if (all(is.finite(final$jacobian)))
    drop(crossprod(final$jacobian, final$moment)) else rep(Inf, parameter_count)
  max_moment <- max(abs(final$moment)); max_score <- max(abs(final_score))
  converged <- is.finite(max_moment) &&
    ((exact && max_moment <= cfg$q_accept_tolerance) ||
     (!exact && is.finite(max_score) && max_score <= cfg$q_accept_tolerance))
  if (!converged) {
    stop(sprintf("Treatment bridge (%s) failed: max moment %.3e, max score %.3e",
                 model, max_moment, max_score), call. = FALSE)
  }
  list(
    model = model, coefficients = parameter,
    basis_scale = b_standard$scale, calibration_scale = g_standard$scale,
    max_abs_moment = max_moment, max_abs_score = max_score,
    jacobian_condition = kappa(final$jacobian, exact = TRUE),
    iterations = iterations, line_search_failures = line_search_failures,
    max_abs_eta_training = max(abs(final$eta)), converged = converged
  )
}

predict_treatment_bridge <- function(fit, data, cfg) {
  basis <- treatment_basis(data, fit$model, proxy = "Z")
  basis <- sweep(basis, 2L, fit$basis_scale, "/")
  eta <- drop(basis %*% fit$coefficients)
  if (any(!is.finite(eta)) || max(eta) >= log(.Machine$double.xmax))
    stop("Treatment-bridge prediction overflow", call. = FALSE)
  q_unclipped <- exp(eta)
  q <- if (is.finite(cfg$q_prediction_cap))
    pmin(q_unclipped, cfg$q_prediction_cap) else q_unclipped
  attr(q, "fraction_clipped") <- if (is.finite(cfg$q_prediction_cap))
    mean(q_unclipped > cfg$q_prediction_cap) else 0
  attr(q, "max_abs_eta") <- max(abs(eta))
  attr(q, "max_q") <- max(q)
  attr(q, "q99") <- unname(quantile(q, 0.99, type = 8))
  attr(q, "ess_fraction") <- (sum(q)^2 / sum(q^2)) / length(q)
  q
}

crossfit_bridges <- function(data, fold, cfg) {
  n <- nrow(data)
  scenarios <- cfg$scenarios
  xi <- matrix(NA_real_, n, nrow(scenarios),
               dimnames = list(NULL, scenarios$scenario))
  h_delta <- matrix(NA_real_, n, 2L,
                    dimnames = list(NULL, c("correct", "omit_w")))
  diagnostics <- vector("list", cfg$folds)

  for (k in seq_len(cfg$folds)) {
    valid <- which(fold == k); train <- which(fold != k)
    training <- data[train, , drop = FALSE]
    validation <- data[valid, , drop = FALSE]
    h_fit <- list(
      correct = fit_outcome_bridge(training, "correct", cfg),
      omit_w = fit_outcome_bridge(training, "omit_w", cfg)
    )
    q_fit <- list(
      quadratic = fit_treatment_bridge(training, "quadratic", cfg),
      main_effects = fit_treatment_bridge(training, "main_effects", cfg)
    )
    h_observed <- lapply(h_fit, predict_outcome_bridge, data = validation)
    q_observed <- lapply(q_fit, predict_treatment_bridge, data = validation, cfg = cfg)

    dose <- outcome_dose_basis(validation$A)
    means <- c(X1 = mean(validation$X1), X2 = mean(validation$X2),
               tanh_W = mean(tanh(validation$W)))
    marginal <- list()
    for (model in names(h_fit)) {
      coefficient <- h_fit[[model]]$coefficients
      value <- coefficient[["(Intercept)"]] + coefficient[["A"]] * dose[, "A"] +
        coefficient[["S_A"]] * dose[, "S_A"] +
        coefficient[["C_A"]] * dose[, "C_A"] +
        coefficient[["X1"]] * means[["X1"]] +
        coefficient[["X2"]] * means[["X2"]]
      delta <- coefficient[["X1"]] * (validation$X1 - means[["X1"]]) +
        coefficient[["X2"]] * (validation$X2 - means[["X2"]])
      if (model == "correct") {
        value <- value + coefficient[["tanh_W"]] * means[["tanh_W"]]
        delta <- delta + coefficient[["tanh_W"]] *
          (tanh(validation$W) - means[["tanh_W"]])
      }
      marginal[[model]] <- value
      h_delta[valid, model] <- delta
    }

    for (row in seq_len(nrow(scenarios))) {
      h_model <- scenarios$outcome_model[[row]]
      q_model <- scenarios$treatment_model[[row]]
      xi[valid, row] <- q_observed[[q_model]] *
        (validation$Y - h_observed[[h_model]]) + marginal[[h_model]]
    }

    diagnostics[[k]] <- data.frame(
      fold = k,
      nuisance = c("h_correct", "h_omit_w", "q_quadratic", "q_main_effects"),
      max_abs_moment = c(h_fit$correct$max_abs_moment, h_fit$omit_w$max_abs_moment,
                         q_fit$quadratic$max_abs_moment,
                         q_fit$main_effects$max_abs_moment),
      max_abs_score = c(h_fit$correct$max_abs_score, h_fit$omit_w$max_abs_score,
                        q_fit$quadratic$max_abs_score,
                        q_fit$main_effects$max_abs_score),
      condition = c(h_fit$correct$condition, h_fit$omit_w$condition,
                    q_fit$quadratic$jacobian_condition,
                    q_fit$main_effects$jacobian_condition),
      iterations = c(NA_integer_, NA_integer_, q_fit$quadratic$iterations,
                     q_fit$main_effects$iterations),
      converged = c(TRUE, TRUE, q_fit$quadratic$converged,
                    q_fit$main_effects$converged),
      fraction_clipped = c(NA_real_, NA_real_,
                           attr(q_observed$quadratic, "fraction_clipped"),
                           attr(q_observed$main_effects, "fraction_clipped")),
      max_abs_eta = c(NA_real_, NA_real_, attr(q_observed$quadratic, "max_abs_eta"),
                      attr(q_observed$main_effects, "max_abs_eta")),
      max_q = c(NA_real_, NA_real_, attr(q_observed$quadratic, "max_q"),
                attr(q_observed$main_effects, "max_q")),
      q99 = c(NA_real_, NA_real_, attr(q_observed$quadratic, "q99"),
              attr(q_observed$main_effects, "q99")),
      ess_fraction = c(NA_real_, NA_real_,
                       attr(q_observed$quadratic, "ess_fraction"),
                       attr(q_observed$main_effects, "ess_fraction")),
      stringsAsFactors = FALSE
    )
  }
  assert_true(all(is.finite(xi)), "Non-finite cross-fitted pseudo-outcome")
  assert_true(all(is.finite(h_delta)), "Non-finite outcome-bridge centering term")
  list(xi = xi, h_delta = h_delta, fold = fold,
       diagnostics = do.call(rbind, diagnostics))
}
