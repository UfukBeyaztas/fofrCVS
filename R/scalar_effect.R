scalar_effect <- function(object,
                          covariate = 1,
                          method = c("local", "CVS", "pCVS"),
                          t_grid = NULL) {
  method <- match.arg(method)
  target_fit <- object$target_fit
  if (target_fit$Q_number < 1) {
    stop("The model has no scalar covariates")
  }
  if (is.character(covariate)) {
    if (length(covariate) != 1 ||
        !covariate %in% target_fit$covariate_names) {
      stop("covariate must name one of: ",
           paste(target_fit$covariate_names, collapse = ", "))
    }
    covariate <- match(covariate, target_fit$covariate_names)
  }
  if (length(covariate) != 1 || !is.finite(covariate) ||
      covariate != as.integer(covariate) || covariate < 1 ||
      covariate > target_fit$Q_number) {
    stop("covariate must be between 1 and ", target_fit$Q_number)
  }
  coefficient <- coefficient_vector(object, method)
  if (is.null(t_grid)) {
    t_grid <- seq(target_fit$y_basis$boundary[1],
                  target_fit$y_basis$boundary[2], length.out = 201)
  }
  Gamma <- matrix(coefficient, target_fit$L, target_fit$My)
  row <- target_fit$P_number * target_fit$Mx + covariate
  Theta_t <- evaluate_basis(target_fit$y_basis, t_grid)
  data.frame(t_grid = t_grid,
             alpha = as.numeric(Theta_t %*% Gamma[row, ]),
             method = method, covariate = covariate)
}
