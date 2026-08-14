coefficient_surface <- function(object,
                                predictor = 1,
                                method = c("local", "CVS", "pCVS"),
                                s_grid = NULL,
                                t_grid = NULL) {
  method <- match.arg(method)
  target_fit <- object$target_fit
  if (is.character(predictor)) {
    if (length(predictor) != 1 ||
        !predictor %in% target_fit$predictor_names) {
      stop("predictor must name one of: ",
           paste(target_fit$predictor_names, collapse = ", "))
    }
    predictor <- match(predictor, target_fit$predictor_names)
  }
  if (length(predictor) != 1 || !is.finite(predictor) ||
      predictor != as.integer(predictor) || predictor < 1 ||
      predictor > target_fit$P_number) {
    stop("predictor must be between 1 and ", target_fit$P_number)
  }
  coefficient <- coefficient_vector(object, method)
  if (is.null(s_grid)) {
    s_grid <- seq(target_fit$x_basis$boundary[1],
                  target_fit$x_basis$boundary[2], length.out = 101)
  }
  if (is.null(t_grid)) {
    t_grid <- seq(target_fit$y_basis$boundary[1],
                  target_fit$y_basis$boundary[2], length.out = 101)
  }
  Gamma <- matrix(coefficient, target_fit$L, target_fit$My)
  rows <- ((predictor - 1) * target_fit$Mx + 1):(predictor * target_fit$Mx)
  Phi_s <- evaluate_basis(target_fit$x_basis, s_grid)
  Theta_t <- evaluate_basis(target_fit$y_basis, t_grid)
  list(s_grid = s_grid, t_grid = t_grid,
       beta = Phi_s %*% Gamma[rows, , drop = FALSE] %*% t(Theta_t),
       method = method, predictor = predictor)
}