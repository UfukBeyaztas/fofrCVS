select_curve_smoothing <- function(Z, Phi, W, rho_grid = NULL) {
  if (is.null(rho_grid)) {
    rho_grid <- default_penalty_grid(crossprod(Phi), W)
  }
  
  n <- nrow(Z)
  J <- ncol(Z)
  scores <- numeric(length(rho_grid))
  
  for (j in seq_along(rho_grid)) {
    P <- solve_psd(crossprod(Phi) + rho_grid[j] * W, t(Phi))
    smoother <- Phi %*% P
    residuals <- Z - Z %*% t(smoother)
    effective_df <- sum(diag(smoother))
    denominator <- max(1 - effective_df / J, 1e-8)
    scores[j] <- mean(residuals^2) / denominator^2
  }
  
  list(
    value = rho_grid[which.min(scores)],
    grid = rho_grid,
    gcv = scores
  )
}