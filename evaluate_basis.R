evaluate_basis <- function(basis, grid) {
  grid <- as.numeric(grid)
  
  if (basis$type == "bspline") {
    return(splineDesign(
      knots = basis$full_knots,
      x = grid,
      ord = basis$degree + 1,
      derivs = rep(0, length(grid)),
      outer.ok = TRUE
    ))
  }
  
  domain_length <- basis$boundary[2] - basis$boundary[1]
  u <- (grid - basis$boundary[1]) / domain_length
  
  Phi <- matrix(0, nrow = length(grid), ncol = basis$M)
  Phi[, 1] <- 1 / sqrt(domain_length)
  
  for (m in 2:basis$M) {
    if (m %% 2 == 0) {
      frequency <- m / 2
      Phi[, m] <- sqrt(2 / domain_length) * cos(2 * pi * frequency * u)
    } else {
      frequency <- (m - 1) / 2
      Phi[, m] <- sqrt(2 / domain_length) * sin(2 * pi * frequency * u)
    }
  }
  
  Phi
}