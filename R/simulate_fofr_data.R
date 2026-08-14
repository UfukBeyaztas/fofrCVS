simulate_fofr_data <- function(n,
                               s_grid = seq(0, 1, length.out = 101),
                               t_grid = seq(0, 1, length.out = 101),
                               covariance_scale = 4,
                               effect = c("common", "shifted"),
                               measurement_sd = 0.05,
                               error_sd = 0.10) {
  effect <- match.arg(effect)
  if (length(n) != 1 || !is.finite(n) || n != as.integer(n) || n < 5) {
    stop("n must be a single integer of at least 5")
  }
  s_grid <- as.numeric(s_grid)
  t_grid <- as.numeric(t_grid)
  if (length(s_grid) < 2 || any(!is.finite(s_grid)) ||
      any(diff(s_grid) <= 0)) {
    stop("s_grid must contain at least two increasing finite values")
  }
  if (length(t_grid) < 2 || any(!is.finite(t_grid)) ||
      any(diff(t_grid) <= 0)) {
    stop("t_grid must contain at least two increasing finite values")
  }
  if (length(covariance_scale) != 1 || !is.finite(covariance_scale) ||
      covariance_scale <= 0) {
    stop("covariance_scale must be a positive finite number")
  }
  if (length(measurement_sd) != 1 || !is.finite(measurement_sd) ||
      measurement_sd < 0 || length(error_sd) != 1 ||
      !is.finite(error_sd) || error_sd < 0) {
    stop("measurement_sd and error_sd must be nonnegative finite numbers")
  }
  
  number_s <- length(s_grid)
  number_t <- length(t_grid)
  covariance_matrix <- covariance_scale *
    exp(-10 * abs(outer(s_grid, s_grid, "-")))
  covariance_factor <- chol(
    covariance_matrix + diag(1e-9, number_s)
  )
  
  latent_activity <- matrix(rnorm(n * number_s), nrow = n) %*%
    covariance_factor
  latent_exposure <- matrix(rnorm(n * number_s), nrow = n) %*%
    covariance_factor
  w <- cbind(
    continuous = rnorm(n),
    binary = rbinom(n, size = 1, prob = 0.5)
  )
  
  beta_activity <- outer(
    s_grid, t_grid,
    function(s, t) sin(2 * pi * s) * cos(pi * t)
  )
  beta_exposure <- outer(
    s_grid, t_grid,
    function(s, t) 2 * (s - 0.5) * (t - 0.5)
  )
  alpha_continuous <- 1 + 0.5 * sin(2 * pi * t_grid)
  alpha_binary <- -0.8 * t_grid
  
  if (effect == "shifted") {
    beta_activity <- -0.35 * beta_activity
    beta_exposure <- beta_exposure + outer(
      s_grid, t_grid,
      function(s, t) 1.25 * cos(pi * s) * sin(2 * pi * t)
    )
    alpha_continuous <- -0.50 * alpha_continuous
    alpha_binary <- alpha_binary + 0.9 * sin(pi * t_grid)
  }
  
  integration_weights <- trapezoid_weights(s_grid)
  conditional_mean <-
    latent_activity %*% (integration_weights * beta_activity) +
    latent_exposure %*% (integration_weights * beta_exposure) +
    w[, "continuous"] %o% alpha_continuous +
    w[, "binary"] %o% alpha_binary
  Y <- conditional_mean +
    matrix(rnorm(n * number_t, sd = error_sd), nrow = n)
  
  observed_activity <- latent_activity +
    matrix(rnorm(n * number_s, sd = measurement_sd), nrow = n)
  observed_exposure <- latent_exposure +
    matrix(rnorm(n * number_s, sd = measurement_sd), nrow = n)
  
  list(
    Y = Y,
    y_grid = t_grid,
    X_list = list(
      activity = observed_activity,
      exposure = observed_exposure
    ),
    x_grids = list(activity = s_grid, exposure = s_grid),
    w = w,
    conditional_mean = conditional_mean,
    latent_X_list = list(
      activity = latent_activity,
      exposure = latent_exposure
    ),
    truth = list(
      beta = list(
        activity = beta_activity,
        exposure = beta_exposure
      ),
      alpha = list(
        continuous = alpha_continuous,
        binary = alpha_binary
      )
    ),
    effect = effect,
    covariance_scale = covariance_scale
  )
}
