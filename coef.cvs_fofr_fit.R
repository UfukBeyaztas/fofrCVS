coef.cvs_fofr_fit <- function(object,
                              method = c("local", "CVS", "pCVS"),
                              ...) {
  method <- match.arg(method)
  target_fit <- object$target_fit
  Gamma <- matrix(
    coefficient_vector(object, method),
    nrow = target_fit$L,
    ncol = target_fit$My
  )
  functional_rows <- unlist(lapply(target_fit$predictor_names, function(name) {
    paste0(name, ".", seq_len(target_fit$Mx))
  }), use.names = FALSE)
  scalar_rows <- target_fit$covariate_names
  rownames(Gamma) <- c(functional_rows, scalar_rows)
  colnames(Gamma) <- paste0("response_basis.", seq_len(target_fit$My))
  attr(Gamma, "method") <- method
  Gamma
}