trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2 || any(diff(grid) <= 0)) {
    stop("grid must contain at least two strictly increasing values")
  }
  
  increments <- diff(grid)
  weights <- numeric(length(grid))
  weights[1] <- increments[1] / 2
  weights[length(grid)] <- increments[length(increments)] / 2
  
  if (length(grid) > 2) {
    weights[2:(length(grid) - 1)] <-
      (increments[1:(length(increments) - 1)] + increments[2:length(increments)]) / 2
  }
  
  weights
}