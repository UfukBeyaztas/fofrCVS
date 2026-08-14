combine_cvs_fofr <- function(target_fit,
                             source_summaries,
                             zeta = NULL,
                             zeta_grid = NULL,
                             target_validation = NULL,
                             tolerance = 1e-7,
                             max_iterations = 20000) {
  if (!inherits(target_fit, "local_fofr_fit")) {
    stop("target_fit must be created by fit_local_fofr")
  }
  if (length(source_summaries) < 1) {
    stop("At least one source summary is required")
  }
  
  summaries <- c(list(fofr_site_summary(target_fit)),
                 lapply(source_summaries, as_fofr_site_summary))
  dims_match <- vapply(summaries, function(x) {
    x$Mx == target_fit$Mx && x$My == target_fit$My &&
      x$P_number == target_fit$P_number && x$Q_number == target_fit$Q_number
  }, logical(1))
  if (!all(dims_match)) {
    stop("All local fits must share the same bases, Mx, My, and predictors")
  }
  signatures_match <- vapply(summaries, function(x) {
    !is.null(x$basis_signature) &&
      isTRUE(all.equal(x$basis_signature, target_fit$basis_signature,
                       tolerance = 1e-12))
  }, logical(1))
  if (!all(signatures_match)) {
    stop("All sites must use identical x/y bases and the same ordered ",
         "functional predictors and scalar covariates")
  }
  if (!is.null(zeta) &&
      (length(zeta) != 1 || !is.finite(zeta) || zeta < 0)) {
    stop("zeta must be a single nonnegative finite number")
  }
  if (!is.null(zeta_grid) &&
      (length(zeta_grid) < 1 || any(!is.finite(zeta_grid)) ||
       any(zeta_grid < 0))) {
    stop("zeta_grid must contain finite nonnegative values")
  }
  if (length(tolerance) != 1 || !is.finite(tolerance) || tolerance <= 0 ||
      length(max_iterations) != 1 || !is.finite(max_iterations) ||
      max_iterations < 1) {
    stop("tolerance and max_iterations must be positive")
  }
  
  group_size <- target_fit$L * target_fit$My
  K <- length(source_summaries)
  
  coefficients <- lapply(summaries, function(x) as.numeric(x$gamma_hat))
  expected_coefficients <- lapply(summaries,
                                  function(x) as.numeric(x$expected_gamma))
  variances <- lapply(summaries, function(x) symmetrize(x$variance_gamma))
  precisions <- lapply(variances, inverse_psd)
  
  pooled_variance <- inverse_psd(Reduce(`+`, precisions))
  source_precision_columns <- do.call(cbind, precisions[-1])
  U_star <- pooled_variance %*% source_precision_columns
  
  delta_hat <- unlist(lapply(seq_len(K), function(k) {
    coefficients[[1]] - coefficients[[k + 1]]
  }), use.names = FALSE)
  expected_delta <- unlist(lapply(seq_len(K), function(k) {
    expected_coefficients[[1]] - expected_coefficients[[k + 1]]
  }), use.names = FALSE)
  
  gamma_cvs <- as.numeric(coefficients[[1]] -
                            U_star %*% (delta_hat - expected_delta))
  
  source_precision_rows <- do.call(rbind, precisions[-1])
  B1 <- block_diagonal(precisions[-1])
  B2 <- source_precision_rows %*% pooled_variance %*% source_precision_columns
  precision_delta <- regularize_psd(symmetrize(B1 - B2))
  
  functional_metric <- coefficient_metric(target_fit)
  metric_root <- matrix_root(functional_metric)
  metric_inverse_root <- matrix_root(functional_metric, inverse = TRUE)
  theta_coefficients <- lapply(coefficients, function(value) {
    as.numeric(metric_root %*% value)
  })
  theta_expected <- lapply(expected_coefficients, function(value) {
    as.numeric(metric_root %*% value)
  })
  theta_variances <- lapply(variances, function(value) {
    regularize_psd(symmetrize(metric_root %*% value %*% metric_root))
  })
  theta_precisions <- lapply(theta_variances, inverse_psd)
  theta_pooled_variance <- inverse_psd(Reduce(`+`, theta_precisions))
  theta_source_precision_columns <- do.call(cbind, theta_precisions[-1])
  theta_U_star <- theta_pooled_variance %*% theta_source_precision_columns
  theta_delta_hat <- unlist(lapply(seq_len(K), function(k) {
    theta_coefficients[[1]] - theta_coefficients[[k + 1]]
  }), use.names = FALSE)
  theta_expected_delta <- unlist(lapply(seq_len(K), function(k) {
    theta_expected[[1]] - theta_expected[[k + 1]]
  }), use.names = FALSE)
  theta_source_precision_rows <- do.call(rbind, theta_precisions[-1])
  theta_B1 <- block_diagonal(theta_precisions[-1])
  theta_B2 <- theta_source_precision_rows %*% theta_pooled_variance %*%
    theta_source_precision_columns
  theta_precision_delta <- regularize_psd(
    symmetrize(theta_B1 - theta_B2)
  )
  
  if (is.null(zeta_grid)) {
    zeta_grid <- penalty_path(
      delta_hat = theta_delta_hat,
      precision_delta = theta_precision_delta,
      group_size = group_size,
      number_groups = K
    )
  }
  
  if (!is.null(zeta)) {
    zeta_grid <- c(zeta, zeta_grid)
  }
  
  zeta_grid <- sort(
    unique(c(zeta_grid, 0)),
    decreasing = TRUE
  )
  
  pcvs_coefficients <- matrix(
    NA_real_,
    nrow = length(zeta_grid),
    ncol = group_size
  )
  
  delta_path <- matrix(
    NA_real_,
    nrow = length(zeta_grid),
    ncol = group_size * K
  )
  
  theta_delta_path <- matrix(
    NA_real_,
    nrow = length(zeta_grid),
    ncol = group_size * K
  )
  convergence <- logical(length(zeta_grid))
  iterations <- integer(length(zeta_grid))
  kkt_residual <- numeric(length(zeta_grid))
  relative_kkt_residual <- numeric(length(zeta_grid))
  current_delta <- rep(0, group_size * K)
  path_lipschitz <- max(
    2 * max(eigen(theta_precision_delta, symmetric = TRUE,
                  only.values = TRUE)$values),
    1e-12
  )
  
  for (j in seq_along(zeta_grid)) {
    solution <- solve_group_lasso(
      delta_hat = theta_delta_hat,
      precision_delta = theta_precision_delta,
      zeta = zeta_grid[j],
      group_size = group_size,
      number_groups = K,
      initial = current_delta,
      tolerance = tolerance,
      max_iterations = max_iterations,
      lipschitz_constant = path_lipschitz
    )
    current_delta <- solution$delta
    theta_delta_path[j, ] <- current_delta
    theta_pcvs <- as.numeric(
      theta_coefficients[[1]] -
        theta_U_star %*% (theta_delta_hat - current_delta)
    )
    pcvs_coefficients[j, ] <- as.numeric(metric_inverse_root %*% theta_pcvs)
    for (k in seq_len(K)) {
      indices <- ((k - 1) * group_size + 1):(k * group_size)
      delta_path[j, indices] <- as.numeric(
        metric_inverse_root %*% current_delta[indices]
      )
    }
    convergence[j] <- solution$converged
    iterations[j] <- solution$iterations
    kkt_residual[j] <- solution$kkt_residual
    relative_kkt_residual[j] <- solution$relative_kkt_residual
  }
  if (any(!convergence)) {
    warning(sum(!convergence), " pCVS path solution(s) did not satisfy the ",
            "convergence and KKT tolerances; these solutions are excluded ",
            "from validation selection")
  }
  
  validation_mse <- rep(NA_real_, length(zeta_grid))
  selected_index <- NULL
  selection <- "not selected"
  
  if (!is.null(target_validation)) {
    validation_data <- validate_fofr_data(target_validation)
    design_validation <- prepare_fofr_design(
      target_fit, validation_data$X_list, validation_data$w,
      validation_data$x_grids
    )
    response_validation <- prepare_response_grid(
      target_fit, validation_data$y_grid
    )
    response_weights <- trapezoid_weights(validation_data$y_grid)
    for (j in seq_along(zeta_grid)) {
      Gamma_j <- matrix(pcvs_coefficients[j, ], target_fit$L, target_fit$My)
      prediction <- sweep(
        design_validation %*% Gamma_j %*%
          t(response_validation$basis_matrix),
        2, response_validation$mean_curve, "+"
      )
      squared_error <- (validation_data$Y - prediction)^2
      validation_mse[j] <- mean(as.numeric(squared_error %*% response_weights))
    }
    validation_mse[!convergence] <- Inf
    if (!any(is.finite(validation_mse))) {
      stop("No pCVS path solution converged; increase max_iterations or ",
           "relax the numerical tolerance")
    }
    selected_index <- which.min(validation_mse)
    selection <- "target validation integrated curve MSE"
  } else if (!is.null(zeta)) {
    selected_index <- which.min(abs(zeta_grid - zeta))
    if (!convergence[selected_index]) {
      stop("The pCVS solution at the requested zeta did not converge")
    }
    selection <- "user-specified zeta"
  }
  
  gamma_pcvs <- if (is.null(selected_index)) NULL else
    as.numeric(pcvs_coefficients[selected_index, ])
  
  structure(list(
    target_fit = target_fit,
    source_summaries = lapply(source_summaries, as_fofr_site_summary),
    gamma_local = coefficients[[1]],
    gamma_cvs = gamma_cvs,
    gamma_pcvs = gamma_pcvs,
    U_star = U_star,
    U_star_metric = theta_U_star,
    delta_hat = delta_hat,
    expected_delta = expected_delta,
    precision_delta = precision_delta,
    pooled_variance = pooled_variance,
    functional_metric = functional_metric,
    metric_root = metric_root,
    metric_inverse_root = metric_inverse_root,
    delta_hat_metric = theta_delta_hat,
    expected_delta_metric = theta_expected_delta,
    precision_delta_metric = theta_precision_delta,
    pooled_variance_metric = theta_pooled_variance,
    zeta_grid = zeta_grid,
    delta_path = delta_path,
    delta_path_metric = theta_delta_path,
    pcvs_coefficient_path = pcvs_coefficients,
    validation_mse = validation_mse,
    selected_zeta = if (is.null(selected_index)) NULL else
      zeta_grid[selected_index],
    selected_index = selected_index,
    selection = selection,
    convergence = convergence,
    iterations = iterations,
    kkt_residual = kkt_residual,
    relative_kkt_residual = relative_kkt_residual,
    K = K,
    group_size = group_size
  ), class = "cvs_fofr_fit")
}
