simulation_config <- function(profile = c("production", "pilot", "test")) {
  profile <- match.arg(profile)
  curve_grid <- seq(-0.9, 0.9, length.out = 61L)

  cfg <- list(
    version = "proximal-joint-level-derivative-v2.0.0",
    profile = profile,
    dgp = "bounded_sinusoidal_nonquadratic",
    sample_sizes = c(2500L, 5000L, 7500L, 10000L),
    repetitions = 5000L,
    base_seed = 20260823L,
    folds = 5L,
    kernel = "epanechnikov",
    curve_grid = curve_grid,
    point_grid = curve_grid[c(1L, 8L, 16L, 24L, 31L, 38L, 46L, 54L, 61L)],
    alpha = 0.05,
    gaussian_band_draws = 10000L,
    level_cv_trim = 0.05,
    cv_min_neighbors = 20L,
    inference_min_neighbors = 8L,
    cv_leverage_limit = 0.95,
    cv1_multipliers = exp(seq(log(0.35), log(2.75), length.out = 15L)),
    cv1_refine = TRUE,
    cv2_j = seq.int(-15L, 15L),
    cv2_log2_step = 1 / 10,
    pi2_multiplier = 1.5,
    level_pi_pilot_multiplier = 1.5,
    derivative_pi_pilot_multiplier = 1.5,
    derivative_pi_pilot_exponent = -1 / 9,
    derivative_rate_exponent = -1 / 7,
    derivative_cv_domain = c(-1, 1),
    derivative_cv_flat_domain = c(-0.9, 0.9),
    derivative_density_floor = 1e-8,
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
    pi_variance_scale = "scale_corrected",
    store_influence_functions = FALSE,
    scenarios = data.frame(
      scenario = c("CC", "CM", "MC", "MM"),
      outcome_model = c("correct", "correct", "omit_w", "omit_w"),
      treatment_model = c("quadratic", "main_effects", "quadratic", "main_effects"),
      stringsAsFactors = FALSE
    ),
    methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2", "ProxLL-US"),
    pdr_methods = c("PDR-PI1", "PDR-PI2", "PDR-CV1", "PDR-CV2"),
    estimands = c("dose_response", "derivative")
  )

  if (profile == "pilot") {
    cfg$sample_sizes <- c(1000L, 2500L)
    cfg$repetitions <- 20L
    cfg$gaussian_band_draws <- 1000L
    cfg$curve_grid <- seq(-0.9, 0.9, length.out = 21L)
    cfg$point_grid <- cfg$curve_grid[c(1L, 4L, 6L, 9L, 11L, 13L, 16L, 18L, 21L)]
    cfg$cv1_multipliers <- exp(seq(log(0.45), log(2.25), length.out = 9L))
  }

  if (profile == "test") {
    cfg$sample_sizes <- 800L
    cfg$repetitions <- 2L
    cfg$folds <- 4L
    cfg$gaussian_band_draws <- 400L
    cfg$curve_grid <- seq(-0.9, 0.9, length.out = 11L)
    cfg$point_grid <- cfg$curve_grid[c(1L, 2L, 3L, 5L, 6L, 7L, 9L, 10L, 11L)]
    cfg$cv1_multipliers <- exp(seq(log(0.5), log(2), length.out = 7L))
    cfg$cv2_j <- seq.int(-8L, 8L)
    cfg$cv_min_neighbors <- 8L
    cfg$inference_min_neighbors <- 6L
  }
  cfg
}
