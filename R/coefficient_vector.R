coefficient_vector <- function(object, method) {
  if (!inherits(object, "cvs_fofr_fit")) {
    stop("object must be created by fit_cvs_fofr or combine_cvs_fofr")
  }
  coefficient <- switch(
    method,
    local = object$gamma_local,
    CVS = object$gamma_cvs,
    pCVS = object$gamma_pcvs
  )
  if (is.null(coefficient)) {
    stop("pCVS has not been selected; provide zeta or target_validation")
  }
  coefficient
}
