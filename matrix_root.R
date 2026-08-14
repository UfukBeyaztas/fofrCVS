matrix_root <- function(A, inverse = FALSE,
                        relative_floor = 1e-10,
                        absolute_floor = 1e-12) {
  A <- symmetrize(as.matrix(A))
  eig <- eigen(A, symmetric = TRUE)
  scale <- max(abs(eig$values))
  floor_value <- max(absolute_floor, relative_floor * scale)
  values <- pmax(eig$values, floor_value)
  weights <- if (inverse) 1 / sqrt(values) else sqrt(values)
  symmetrize(eig$vectors %*% (t(eig$vectors) * weights))
}