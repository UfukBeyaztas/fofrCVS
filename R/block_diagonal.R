block_diagonal <- function(matrices) {
  dimensions <- vapply(matrices, nrow, numeric(1))
  answer <- matrix(0, sum(dimensions), sum(dimensions))
  starts <- cumsum(c(1, dimensions[-length(dimensions)]))
  
  for (j in seq_along(matrices)) {
    indices <- starts[j]:(starts[j] + dimensions[j] - 1)
    answer[indices, indices] <- matrices[[j]]
  }
  
  answer
}
