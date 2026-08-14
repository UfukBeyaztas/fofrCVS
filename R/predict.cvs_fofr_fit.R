predict.cvs_fofr_fit <- function(object,
                                 newdata,
                                 method = c("local", "CVS", "pCVS"),
                                 ...) {
  method <- match.arg(method)
  coefficient <- coefficient_vector(object, method)
  if (!is.list(newdata) || is.null(newdata$X_list)) {
    stop("newdata must be a list containing X_list, x_grids, and optionally w")
  }
  target_fit <- object$target_fit
  design <- prepare_fofr_design(
    target_fit, newdata$X_list, newdata$w, newdata$x_grids
  )
  response_grid <- prepare_response_grid(target_fit, newdata$y_grid)
  Gamma <- matrix(coefficient, target_fit$L, target_fit$My)
  sweep(design %*% Gamma %*% t(response_grid$basis_matrix),
        2, response_grid$mean_curve, "+")
}
