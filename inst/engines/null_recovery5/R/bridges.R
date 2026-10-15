quadratic_basis <- function(proxy, a, x1, x2) {
  cbind(
    `(Intercept)` = 1, proxy = proxy, A = a, X1 = x1, X2 = x2,
    proxy2 = proxy^2, A2 = a^2, X1_2 = x1^2, X2_2 = x2^2,
    proxy_A = proxy * a, proxy_X1 = proxy * x1, proxy_X2 = proxy * x2,
    A_X1 = a * x1, A_X2 = a * x2, X1_X2 = x1 * x2
  )
}

main_effects_basis <- function(proxy, a, x1, x2) {
  cbind(`(Intercept)` = 1, proxy = proxy, A = a, X1 = x1, X2 = x2)
}

treatment_basis <- function(data, model = c("quadratic", "main_effects"),
                            proxy = c("Z", "W")) {
  model <- match.arg(model); proxy <- match.arg(proxy)
  if (model == "quadratic")
    quadratic_basis(data[[proxy]], data$A, data$X1, data$X2)
  else main_effects_basis(data[[proxy]], data$A, data$X1, data$X2)
}

product_calibration_target <- function(data) {
  ma <- mean(data$A); mx1 <- mean(data$X1); mx2 <- mean(data$X2); mw <- mean(data$W)
  c(
    `(Intercept)` = 1, proxy = mw, A = ma, X1 = mx1, X2 = mx2,
    proxy2 = mean(data$W^2), A2 = mean(data$A^2),
    X1_2 = mean(data$X1^2), X2_2 = mean(data$X2^2),
    proxy_A = mw * ma, proxy_X1 = mean(data$W * data$X1),
    proxy_X2 = mean(data$W * data$X2), A_X1 = ma * mx1,
    A_X2 = ma * mx2, X1_X2 = mean(data$X1 * data$X2)
  )
}

outcome_dose_basis <- function(a) {
  variance_a <- exposure_variance()
  cbind(`(Intercept)` = 1, A = a, A2c = a^2 - variance_a,
        A3c = a^3 - 3 * variance_a * a)
}

outcome_design <- function(data, model = c("correct", "omit_w"),
                           a = data$A, w = data$W,
                           x1 = data$X1, x2 = data$X2) {
  model <- match.arg(model)
  base <- cbind(outcome_dose_basis(a), X1 = x1, X2 = x2)
  if (model == "correct") cbind(base, tanh_W = tanh(w)) else base
}

outcome_instruments <- function(data) {
  cbind(outcome_dose_basis(data$A), X1 = data$X1, X2 = data$X2, Z = data$Z)
}

fit_outcome_bridge_matrix <- function(data, y, model = c("correct", "omit_w"), cfg) {
  model <- match.arg(model); y <- as.matrix(y)
  hmat <- outcome_design(data, model); zmat <- outcome_instruments(data)
  omega_inverse <- safe_solve(crossprod(zmat) / nrow(data),
                              ridge = cfg$nuisance_ridge,
                              condition_limit = cfg$nuisance_condition_limit)
  g <- crossprod(zmat, hmat) / nrow(data)
  gy <- crossprod(zmat, y) / nrow(data)
  lhs <- crossprod(g, omega_inverse %*% g)
  coefficient <- safe_solve(lhs, crossprod(g, omega_inverse %*% gy),
                            ridge = cfg$nuisance_ridge,
                            condition_limit = cfg$nuisance_condition_limit)
  rownames(coefficient) <- colnames(hmat); colnames(coefficient) <- colnames(y)
  moment <- crossprod(zmat, y - hmat %*% coefficient) / nrow(data)
  score <- crossprod(g, omega_inverse %*% moment)
  list(model = model, coefficients = coefficient,
       max_abs_moment = apply(abs(moment), 2L, max),
       max_abs_score = apply(abs(score), 2L, max),
       condition = kappa(lhs, exact = TRUE))
}

predict_outcome_bridge_matrix <- function(fit, data, a = data$A) {
  outcome_design(data, fit$model, a = a) %*% fit$coefficients
}

standardize_columns <- function(x) {
  scale <- sqrt(colMeans(x^2)); scale[!is.finite(scale) | scale < 1e-8] <- 1
  scale[[1L]] <- 1; list(x = sweep(x, 2L, scale, "/"), scale = scale)
}

fit_treatment_bridge <- function(data, model = c("quadratic", "main_effects"), cfg) {
  model <- match.arg(model)
  b_raw <- treatment_basis(data, model, proxy = "Z")
  g_raw <- treatment_basis(data, "quadratic", proxy = "W")
  target_raw <- product_calibration_target(data)
  assert_true(identical(colnames(g_raw), names(target_raw)), "Calibration basis mismatch")
  b_standard <- standardize_columns(b_raw); g_standard <- standardize_columns(g_raw)
  b <- b_standard$x; g <- g_standard$x; target <- target_raw / g_standard$scale
  parameter_count <- ncol(b); moment_count <- ncol(g)
  evaluate <- function(parameter, jacobian = TRUE) {
    eta <- drop(b %*% parameter)
    if (any(!is.finite(eta)) || max(abs(eta)) > cfg$q_eta_limit)
      return(list(moment = rep(Inf, moment_count), jacobian = matrix(Inf, moment_count,
                                                                     parameter_count)))
    q <- exp(eta); moment <- colMeans(g * q) - target
    list(moment = moment,
         jacobian = if (jacobian) crossprod(g, b * q) / nrow(data) else NULL)
  }
  parameter <- rep(0, parameter_count); iterations <- 0L; line_search_failures <- 0L
  for (iteration in seq_len(cfg$q_max_iter)) {
    iterations <- iteration; current <- evaluate(parameter); norm0 <- sqrt(sum(current$moment^2))
    score0 <- if (all(is.finite(current$jacobian)))
      drop(crossprod(current$jacobian, current$moment)) else rep(Inf, parameter_count)
    exact <- moment_count == parameter_count
    if (is.finite(norm0) && ((exact && max(abs(current$moment)) <= cfg$q_root_tolerance) ||
        (!exact && max(abs(score0)) <= cfg$q_root_tolerance))) break
    step <- tryCatch(if (exact) {
      drop(qr.solve(current$jacobian, current$moment, tol = 1e-12))
    } else {
      normal <- crossprod(current$jacobian)
      drop(qr.solve(normal + diag(cfg$nuisance_ridge, parameter_count),
                    crossprod(current$jacobian, current$moment), tol = 1e-12))
    }, error = function(e) rep(NA_real_, parameter_count))
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
  final <- evaluate(parameter); final_score <- drop(crossprod(final$jacobian, final$moment))
  exact <- moment_count == parameter_count
  needs_fallback <- !all(is.finite(final$moment)) ||
    (exact && max(abs(final$moment)) > cfg$q_accept_tolerance) ||
    (!exact && max(abs(final_score)) > cfg$q_accept_tolerance)
  if (needs_fallback) {
    objective <- function(value) { z <- evaluate(value); if (any(!is.finite(z$moment))) 1e100 else sum(z$moment^2) }
    gradient <- function(value) { z <- evaluate(value); if (any(!is.finite(z$moment))) rep(0, parameter_count) else drop(2 * crossprod(z$jacobian, z$moment)) }
    optimization <- optim(parameter, objective, gradient, method = "BFGS",
                          control = list(maxit = 500L, reltol = 1e-12))
    if (is.finite(optimization$value) && optimization$value < objective(parameter))
      parameter <- optimization$par
    final <- evaluate(parameter); final_score <- drop(crossprod(final$jacobian, final$moment))
  }
  max_moment <- max(abs(final$moment)); max_score <- max(abs(final_score))
  converged <- is.finite(max_moment) && ((exact && max_moment <= cfg$q_accept_tolerance) ||
    (!exact && is.finite(max_score) && max_score <= cfg$q_accept_tolerance))
  if (!converged) stop(sprintf("Treatment bridge (%s) failed: moment %.3e, score %.3e",
                                model, max_moment, max_score), call. = FALSE)
  list(model = model, coefficients = parameter, basis_scale = b_standard$scale,
       calibration_scale = g_standard$scale, max_abs_moment = max_moment,
       max_abs_score = max_score, jacobian_condition = kappa(final$jacobian, exact = TRUE),
       iterations = iterations, line_search_failures = line_search_failures,
       converged = converged)
}

predict_treatment_bridge <- function(fit, data, cfg) {
  b <- sweep(treatment_basis(data, fit$model, proxy = "Z"), 2L, fit$basis_scale, "/")
  eta <- drop(b %*% fit$coefficients)
  if (any(!is.finite(eta)) || max(eta) >= log(.Machine$double.xmax))
    stop("Treatment-bridge prediction overflow", call. = FALSE)
  q_raw <- exp(eta); q <- if (is.finite(cfg$q_prediction_cap))
    pmin(q_raw, cfg$q_prediction_cap) else q_raw
  attr(q, "fraction_clipped") <- if (is.finite(cfg$q_prediction_cap))
    mean(q_raw > cfg$q_prediction_cap) else 0
  attr(q, "max_abs_eta") <- max(abs(eta)); attr(q, "max_q") <- max(q)
  attr(q, "q99") <- unname(quantile(q, .99, type = 8))
  attr(q, "ess_fraction") <- (sum(q)^2 / sum(q^2)) / length(q); q
}

crossfit_all_conditions <- function(base, outcome, fold, cfg) {
  n <- nrow(base); conditions <- outcome$conditions; scenarios <- cfg$scenarios
  key <- expand.grid(condition = conditions$condition, scenario = scenarios$scenario,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  xi <- matrix(NA_real_, n, nrow(key), dimnames = list(NULL, paste(key$condition, key$scenario, sep = "__")))
  # For the present additive bridge class,
  # h(W_r,a,X_r)-P_{n,k}^{WX}h(W,a,X) is independent of the outer dose a.
  # Storing this exact same-fold row factor provides a response-level interface
  # for the product-measure Hajek perturbation without materializing an
  # n_k-by-n_k cross-prediction matrix.
  row_calibration_factor <- matrix(
    NA_real_, n, nrow(key), dimnames = dimnames(xi))
  diagnostics <- vector("list", cfg$folds)
  for (k in seq_len(cfg$folds)) {
    valid <- which(fold == k); train <- which(fold != k)
    training <- base[train, , drop = FALSE]; validation <- base[valid, , drop = FALSE]
    y_train <- outcome$y[train, , drop = FALSE]; y_valid <- outcome$y[valid, , drop = FALSE]
    h_fit <- list(correct = fit_outcome_bridge_matrix(training, y_train, "correct", cfg),
                  omit_w = fit_outcome_bridge_matrix(training, y_train, "omit_w", cfg))
    q_fit <- list(quadratic = fit_treatment_bridge(training, "quadratic", cfg),
                  main_effects = fit_treatment_bridge(training, "main_effects", cfg))
    h_observed <- lapply(h_fit, predict_outcome_bridge_matrix, data = validation)
    q_observed <- lapply(q_fit, predict_treatment_bridge, data = validation, cfg = cfg)
    dose <- outcome_dose_basis(validation$A)
    means <- c(X1 = mean(validation$X1), X2 = mean(validation$X2),
               tanh_W = mean(tanh(validation$W)))
    marginal <- list()
    for (model in names(h_fit)) {
      coefficient <- h_fit[[model]]$coefficients
      value <- dose %*% coefficient[colnames(dose), , drop = FALSE]
      offset <- means[["X1"]] * coefficient["X1", ] +
        means[["X2"]] * coefficient["X2", ]
      if (model == "correct") offset <- offset +
        means[["tanh_W"]] * coefficient["tanh_W", ]
      value <- sweep(value, 2L, offset, "+")
      marginal[[model]] <- value
    }
    for (s in seq_len(nrow(scenarios))) {
      hm <- scenarios$outcome_model[[s]]; qm <- scenarios$treatment_model[[s]]
      value <- q_observed[[qm]] * (y_valid - h_observed[[hm]]) + marginal[[hm]]
      row_effect <- h_observed[[hm]] - marginal[[hm]]
      columns <- which(key$scenario == scenarios$scenario[[s]])
      match_condition <- match(key$condition[columns], conditions$condition)
      xi[valid, columns] <- value[, match_condition, drop = FALSE]
      row_calibration_factor[valid, columns] <-
        row_effect[, match_condition, drop = FALSE]
    }
    q_tail <- function(model, field) unname(attr(q_observed[[model]], field))
    diagnostics[[k]] <- data.frame(
      fold = k, nuisance = c("h_correct", "h_omit_w", "q_quadratic", "q_main_effects"),
      fitting_criterion = c("GMM score", "GMM score", "moment root", "GMM score"),
      max_abs_moment = c(max(h_fit$correct$max_abs_moment), max(h_fit$omit_w$max_abs_moment),
                         q_fit$quadratic$max_abs_moment, q_fit$main_effects$max_abs_moment),
      max_abs_score = c(max(h_fit$correct$max_abs_score), max(h_fit$omit_w$max_abs_score),
                        q_fit$quadratic$max_abs_score, q_fit$main_effects$max_abs_score),
      fit_condition = c(h_fit$correct$condition, h_fit$omit_w$condition,
                        q_fit$quadratic$jacobian_condition,
                        q_fit$main_effects$jacobian_condition),
      iterations = c(NA_integer_, NA_integer_, q_fit$quadratic$iterations,
                     q_fit$main_effects$iterations),
      line_search_failures = c(NA_integer_, NA_integer_,
                               q_fit$quadratic$line_search_failures,
                               q_fit$main_effects$line_search_failures),
      max_abs_eta = c(NA, NA, q_tail("quadratic", "max_abs_eta"),
                      q_tail("main_effects", "max_abs_eta")),
      q99 = c(NA, NA, q_tail("quadratic", "q99"),
              q_tail("main_effects", "q99")),
      max_q = c(NA, NA, q_tail("quadratic", "max_q"),
                q_tail("main_effects", "max_q")),
      ess_fraction = c(NA, NA, q_tail("quadratic", "ess_fraction"),
                       q_tail("main_effects", "ess_fraction")),
      fraction_clipped = c(NA, NA, q_tail("quadratic", "fraction_clipped"),
                           q_tail("main_effects", "fraction_clipped")),
      converged = c(TRUE, TRUE, q_fit$quadratic$converged,
                    q_fit$main_effects$converged),
      stringsAsFactors = FALSE)
  }
  assert_true(all(is.finite(xi)), "Non-finite cross-fitted pseudo-outcome")
  assert_true(all(is.finite(row_calibration_factor)),
              "Non-finite row-calibration factor")
  fold_mean_error <- max(vapply(seq_len(cfg$folds), function(k)
    max(abs(colMeans(row_calibration_factor[fold == k, , drop = FALSE]))),
    numeric(1L)))
  assert_true(fold_mean_error <= 2e-10,
              sprintf("Row-Hajek factors are not fold centered (%.3e)",
                      fold_mean_error))
  list(xi = xi, key = key, conditions = conditions, fold = fold,
       row_calibration_factor = row_calibration_factor,
       row_calibration_interface = "low_rank_exact",
       diagnostics = do.call(rbind, diagnostics))
}
