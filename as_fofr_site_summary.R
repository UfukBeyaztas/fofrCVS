as_fofr_site_summary <- function(object) {
  if (inherits(object, "local_fofr_fit")) {
    return(fofr_site_summary(object))
  }
  if (!inherits(object, "cvs_fofr_summary")) {
    stop("Sources must be local_fofr_fit or cvs_fofr_summary objects")
  }
  object
}