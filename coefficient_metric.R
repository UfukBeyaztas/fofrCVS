coefficient_metric <- function(local_fit) {
  row_blocks <- c(
    rep(list(local_fit$x_basis$Psi), local_fit$P_number),
    if (local_fit$Q_number > 0) list(diag(local_fit$Q_number)) else NULL
  )
  row_metric <- block_diagonal(row_blocks)
  regularize_psd(symmetrize(kronecker(local_fit$y_basis$Psi, row_metric)))
}