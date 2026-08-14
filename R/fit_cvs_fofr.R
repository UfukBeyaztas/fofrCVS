fit_cvs_fofr <- function(datasets,
                         x_basis_type = c("fourier", "bspline"),
                         y_basis_type = c("fourier", "bspline"),
                         Mx = 7,
                         My = 5,
                         degree = 3,
                         x_boundary = c(0, 1),
                         y_boundary = c(0, 1),
                         rho_x = NULL,
                         rho_y = NULL,
                         lambda = NULL,
                         lambda_grid = NULL,
                         variance_method = c("HC", "model"),
                         zeta = NULL,
                         zeta_grid = NULL,
                         target_validation = NULL,
                         tolerance = 1e-7,
                         max_iterations = 20000) {
  x_basis_type <- match.arg(x_basis_type)
  y_basis_type <- match.arg(y_basis_type)
  variance_method <- match.arg(variance_method)
  
  if (!is.list(datasets) || length(datasets) < 2) {
    stop("datasets must contain one target followed by at least one source")
  }
  
  x_basis <- fofr_basis(x_basis_type, M = Mx, degree = degree,
                        boundary = x_boundary)
  y_basis <- fofr_basis(y_basis_type, M = My, degree = degree,
                        boundary = y_boundary)
  
  local_fits <- lapply(datasets, function(dataset) {
    fit_local_fofr(dataset, x_basis, y_basis,
                   rho_x = rho_x, rho_y = rho_y, lambda = lambda,
                   lambda_grid = lambda_grid,
                   variance_method = variance_method)
  })
  
  answer <- combine_cvs_fofr(
    target_fit = local_fits[[1]],
    source_summaries = lapply(local_fits[-1], fofr_site_summary),
    zeta = zeta,
    zeta_grid = zeta_grid,
    target_validation = target_validation,
    tolerance = tolerance,
    max_iterations = max_iterations
  )
  answer$local_fits <- local_fits
  answer$call <- match.call()
  answer
}
