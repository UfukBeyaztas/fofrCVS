prepare_response_grid <- function(target_fit, y_grid = NULL) {
  if (is.null(y_grid)) {
    y_grid <- target_fit$y_grid
  }
  y_grid <- as.numeric(y_grid)
  if (length(y_grid) < 2 || any(!is.finite(y_grid)) ||
      any(diff(y_grid) <= 0)) {
    stop("The response grid must be finite and strictly increasing")
  }
  training_grid <- target_fit$y_grid
  grid_tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(training_grid)))
  if (min(y_grid) < min(training_grid) - grid_tolerance ||
      max(y_grid) > max(training_grid) + grid_tolerance) {
    stop("The response grid is outside the target training-grid range")
  }
  same_grid <- length(y_grid) == length(training_grid) &&
    max(abs(y_grid - training_grid)) <= grid_tolerance
  mean_curve <- if (same_grid) {
    target_fit$y_mean_curve
  } else {
    approx(training_grid, target_fit$y_mean_curve,
           xout = y_grid, rule = 2)$y
  }
  list(
    grid = y_grid,
    basis_matrix = evaluate_basis(target_fit$y_basis, y_grid),
    mean_curve = mean_curve
  )
}
