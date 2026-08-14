solve_psd <- function(A, b, relative_floor = 1e-10) {
  inverse_psd(A, relative_floor = relative_floor) %*% b
}
