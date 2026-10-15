gauss_legendre_rule <- function(order, lower, upper) {
  stopifnot(order >= 2L, is.finite(lower), is.finite(upper), lower < upper)
  j <- seq_len(order - 1L)
  off <- j / sqrt(4 * j^2 - 1)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(j, j + 1L)] <- off
  jacobi[cbind(j + 1L, j)] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  node <- eig$values[ord]
  weight <- 2 * eig$vectors[1L, ord]^2
  list(
    node = (upper + lower) / 2 + (upper - lower) * node / 2,
    weight = (upper - lower) * weight / 2
  )
}

composite_gauss_legendre_rule <- function(panels, nodes_per_panel,
                                          lower, upper) {
  stopifnot(panels >= 1L, nodes_per_panel >= 2L, lower < upper)
  boundary <- seq(lower, upper, length.out = panels + 1L)
  pieces <- lapply(seq_len(panels), function(j)
    gauss_legendre_rule(nodes_per_panel, boundary[[j]], boundary[[j + 1L]]))
  list(
    node = unlist(lapply(pieces, `[[`, "node"), use.names = FALSE),
    weight = unlist(lapply(pieces, `[[`, "weight"), use.names = FALSE)
  )
}

simulation_config <- function(profile = c("production", "calibration", "test")) {
  profile <- match.arg(profile)
  # The test integral is distinct from the unchanged 61-point grid used by
  # PI1's pilot-integral calculation and by displayed curves. Composite
  # quadrature resolves the many kernel-induced knots more reliably than a
  # single global rule of the same size.
  quadrature_panels <- 64L
  quadrature_nodes_per_panel <- 8L
  quadrature_order <- quadrature_panels * quadrature_nodes_per_panel
  quadrature <- composite_gauss_legendre_rule(
    quadrature_panels, quadrature_nodes_per_panel, -0.9, 0.9)
  test_weight <- (1 - (quadrature$node / 0.9)^2)^2

  cfg <- list(
    version = "failed5-local-final-cv-admissibility-v1.0.0",
    fold_cv_filter = TRUE,
    fold_cv_filter_rule = "original_final_candidate_cloud_all_fold_actual_operators_raw_and_US",
    pi1_pilot_rule = "original_unfiltered_LL_CV_pilot_unchanged",
    profile = profile,
    sample_sizes = c(2500L, 5000L, 7500L, 10000L),
    repetitions = 1000L,
    bootstrap_draws = 4999L,
    bootstrap_block = 100L,
    bootstrap_apply_backend = "blas",
    multiplier_laws = "two_point",
    bootstrap_schemes = c("fold_parametric", "fold_nonparametric"),
    nonparametric_residual_fit = "same_method_in_sample_same_bandwidth",
    singular_observation_fit = "zero_combined_operator",
    base_seed = 20260825L,
    folds = 5L,
    alpha = 0.05,
    kernel = "epanechnikov",
    null_degrees = 0:2,
    null_labels = c(`0` = "constant", `1` = "linear", `2` = "quadratic"),
    # Departure grid specified by the user on 2026-09-26.
    delta_grid = c(0, 0.004, 0.008, 0.012, 0.016, 0.024, 0.032, 0.040,
                   0.052, 0.064, 0.080),
    delta_grid_status = "predeclared_common_cubic",
    delta_calibration = list(
      config_id = "not_applicable_common_cubic", seed = NA_integer_,
      repetitions = 0L, bootstrap_draws = 0L,
      decision = "USER_SPECIFIED_COMMON_CUBIC_GRID_20260926",
      decision_rule = paste(
        "Use the eleven departure values specified by the user on 2026-09-26.",
        "Keep the DGP equations, nuisance solvers, bandwidth selection,",
        "sample sizes, original seed map and bootstrap settings unchanged."
      )
    ),
    quadrature_order = quadrature_order,
    quadrature_panels = quadrature_panels,
    quadrature_nodes_per_panel = quadrature_nodes_per_panel,
    quadrature_node = quadrature$node,
    quadrature_weight = quadrature$weight,
    test_weight = test_weight,
    weighted_quadrature = quadrature$weight * test_weight,
    curve_grid = seq(-0.9, 0.9, length.out = 61L),
    test_support = c(-0.9, 0.9),
    exposure_variance = 0.8203125,
    methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                "ProxLL-raw", "ProxLL-raw-matched",
                "ProxLL-US", "ProxLL-US-matched"),
    bandwidth_methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2",
                          "ProxLL-raw", "ProxLL-US"),
    pdr_methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2"),
    scenarios = data.frame(
      scenario = c("CC", "CM", "MC", "MM"),
      outcome_model = c("correct", "correct", "omit_w", "omit_w"),
      treatment_model = c("quadratic", "main_effects", "quadratic", "main_effects"),
      identified = c(TRUE, TRUE, TRUE, FALSE),
      stringsAsFactors = FALSE
    ),
    level_cv_trim = 0.05,
    cv_min_neighbors = 20L,
    inference_min_neighbors = 8L,
    fold_inference_min_neighbors = 3L,
    cv_leverage_limit = 0.95,
    cv1_multipliers = exp(seq(log(0.35), log(2.75), length.out = 15L)),
    cv1_refine = TRUE,
    cv2_j = seq.int(-15L, 15L),
    cv2_log2_step = 1 / 10,
    pi2_multiplier = 1.5,
    level_pi_pilot_multiplier = 1.5,
    variance_knn_power = 0.5,
    q_prediction_cap = Inf,
    q_eta_limit = 25,
    q_max_iter = 80L,
    q_root_tolerance = 1e-8,
    q_accept_tolerance = 5e-5,
    nuisance_ridge = 1e-10,
    nuisance_condition_limit = 1e12,
    local_condition_limit = 1e12,
    local_reproduction_tolerance = 2e-8,
    store_bootstrap_statistics = FALSE
  )

  if (profile == "calibration") {
    # The calibration sample is independent of every production replication;
    # otherwise data-adaptive selection of the common delta grid would reuse
    # the first production Monte Carlo samples.
    cfg$base_seed <- 731928641L
    cfg$sample_sizes <- c(2500L, 7500L)
    cfg$repetitions <- 100L
    cfg$bootstrap_draws <- 399L
    cfg$bootstrap_block <- 100L
    cfg$scenarios <- cfg$scenarios[cfg$scenarios$identified, , drop = FALSE]
  }
  if (profile == "test") {
    cfg$base_seed <- 14921731L
    cfg$sample_sizes <- c(800L, 1200L)
    cfg$repetitions <- 2L
    cfg$bootstrap_draws <- 39L
    cfg$bootstrap_block <- 20L
    cfg$delta_grid <- c(0, 0.04)
    cfg$cv1_multipliers <- exp(seq(log(0.5), log(2), length.out = 7L))
    cfg$cv2_j <- seq.int(-8L, 8L)
    cfg$cv_min_neighbors <- 8L
    cfg$inference_min_neighbors <- 6L
  }
  if (profile == "production" &&
      cfg$delta_grid_status != "predeclared_common_cubic") {
    stop("Production is locked until the common-cubic grid is predeclared.",
         call. = FALSE)
  }
  cfg
}
