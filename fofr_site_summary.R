fofr_site_summary <- function(local_fit) {
  if (!inherits(local_fit, "local_fofr_fit")) {
    stop("local_fit must be created by fit_local_fofr")
  }
  structure(list(
    gamma_hat = local_fit$gamma_hat,
    expected_gamma = local_fit$expected_gamma,
    variance_gamma = local_fit$variance_gamma,
    Mx = local_fit$Mx,
    My = local_fit$My,
    P_number = local_fit$P_number,
    Q_number = local_fit$Q_number,
    basis_signature = local_fit$basis_signature
  ), class = "cvs_fofr_summary")
}