print.cvs_fofr_fit <- function(x, ...) {
  cat("Control-variates transfer learning, function-on-function regression\n")
  cat("Target sample size:", x$target_fit$n, "\n")
  cat("Functional predictors:", x$target_fit$P_number,
      " Scalar covariates:", x$target_fit$Q_number, "\n")
  cat("Mx:", x$target_fit$Mx, " My:", x$target_fit$My,
      " Number of sources:", x$K, "\n")
  cat("Target lambda:", format(x$target_fit$lambda, digits = 5), "\n")
  cat("Target design rank:", x$target_fit$design_rank, "of",
      x$target_fit$L, " Condition number:",
      format(x$target_fit$design_condition, digits = 5), "\n")
  if (is.null(x$selected_zeta)) {
    cat("pCVS: no zeta selected\n")
  } else {
    cat("Selected pCVS zeta:", format(x$selected_zeta, digits = 5), "\n")
    cat("Selection rule:", x$selection, "\n")
    cat("pCVS metric: integrated L2 norm of coefficient functions\n")
    cat("Selected solution relative KKT residual:",
        format(x$relative_kkt_residual[x$selected_index], digits = 4), "\n")
  }
  invisible(x)
}