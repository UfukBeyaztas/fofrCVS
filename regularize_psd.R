regularize_psd <- function(A, relative_floor = 1e-8, absolute_floor = 1e-12) {
  A <- symmetrize(as.matrix(A))
  eig <- eigen(A, symmetric = TRUE)
  spectral_scale <- max(abs(eig$values))
  floor_value <- max(absolute_floor, relative_floor * spectral_scale)
  values <- pmax(eig$values, floor_value)
  symmetrize(eig$vectors %*% (t(eig$vectors) * values))
}