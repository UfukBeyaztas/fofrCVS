fit_local_fofr <- function(dataset,
                           x_basis,
                           y_basis,
                           rho_x = NULL,
                           rho_y = NULL,
                           lambda = NULL,
                           lambda_grid = NULL,
                           variance_method = c("HC", "model")) {
  variance_method <- match.arg(variance_method)
  if (!inherits(x_basis, "fofr_basis") ||
      !inherits(y_basis, "fofr_basis")) {
    stop("x_basis and y_basis must be created by fofr_basis")
  }
  dataset <- validate_fofr_data(dataset)
  
  validate_penalty <- function(value, name, allow_vector = FALSE) {
    if (is.null(value)) return(invisible(NULL))
    if ((!allow_vector && length(value) != 1) || any(!is.finite(value)) ||
        any(value < 0)) {
      stop(name, " must contain finite nonnegative values")
    }
  }
  validate_penalty(rho_x, "rho_x", allow_vector = TRUE)
  validate_penalty(rho_y, "rho_y")
  validate_penalty(lambda, "lambda")
  if (!is.null(lambda_grid) &&
      (length(lambda_grid) < 1 || any(!is.finite(lambda_grid)) ||
       any(lambda_grid < 0))) {
    stop("lambda_grid must contain finite nonnegative values")
  }
  
  Y <- dataset$Y
  y_grid <- dataset$y_grid
  X_list <- dataset$X_list
  x_grids <- dataset$x_grids
  w <- dataset$w
  predictor_names <- dataset$predictor_names
  covariate_names <- dataset$covariate_names
  
  n <- nrow(Y)
  P_number <- length(X_list)
  Q_number <- if (is.null(w)) 0 else ncol(w)
  Mx <- x_basis$M
  My <- y_basis$M
  L <- P_number * Mx + Q_number
  
  boundary_tolerance_x <- 100 * .Machine$double.eps *
    max(1, max(abs(x_basis$boundary)))
  boundary_tolerance_y <- 100 * .Machine$double.eps *
    max(1, max(abs(y_basis$boundary)))
  if (min(y_grid) < y_basis$boundary[1] - boundary_tolerance_y ||
      max(y_grid) > y_basis$boundary[2] + boundary_tolerance_y) {
    stop("y_grid must lie inside the response-basis boundary")
  }
  for (p in seq_len(P_number)) {
    if (min(x_grids[[p]]) < x_basis$boundary[1] - boundary_tolerance_x ||
        max(x_grids[[p]]) > x_basis$boundary[2] + boundary_tolerance_x) {
      stop("Functional predictor ", p,
           " has a grid outside the predictor-basis boundary")
    }
  }
  
  y_mean_curve <- colMeans(Y)
  Y_centered <- sweep(Y, 2, y_mean_curve, "-")
  x_mean_curves <- lapply(X_list, colMeans)
  X_centered <- lapply(seq_len(P_number), function(p) {
    sweep(X_list[[p]], 2, x_mean_curves[[p]], "-")
  })
  w_means <- if (Q_number > 0) colMeans(w) else numeric(0)
  w_centered <- if (Q_number > 0) sweep(w, 2, w_means, "-") else NULL
  
  expand_rho <- function(value, count, name) {
    if (is.null(value)) return(rep(list(NULL), count))
    if (length(value) == 1) return(as.list(rep(value, count)))
    if (length(value) != count) {
      stop(name, " must be NULL, a scalar, or have one value per predictor")
    }
    as.list(value)
  }
  rho_x_values <- expand_rho(rho_x, P_number, "rho_x")
  
  Phi_x <- vector("list", P_number)
  P_x <- vector("list", P_number)
  rho_x_used <- numeric(P_number)
  measurement_variances <- numeric(P_number)
  for (p in seq_len(P_number)) {
    Phi_x[[p]] <- evaluate_basis(x_basis, x_grids[[p]])
    if (is.null(rho_x_values[[p]])) {
      sel <- select_curve_smoothing(X_centered[[p]], Phi_x[[p]], x_basis$W)
      rho_x_used[p] <- sel$value
    } else {
      rho_x_used[p] <- rho_x_values[[p]]
    }
    P_x[[p]] <- solve_psd(crossprod(Phi_x[[p]]) + rho_x_used[p] * x_basis$W,
                          t(Phi_x[[p]]))
    smoother_x <- Phi_x[[p]] %*% P_x[[p]]
    fitted_x <- X_centered[[p]] %*% t(smoother_x)
    measurement_variances[p] <- mean((X_centered[[p]] - fitted_x)^2)
  }
  
  Theta <- evaluate_basis(y_basis, y_grid)
  if (is.null(rho_y)) {
    sel_y <- select_curve_smoothing(Y_centered, Theta, y_basis$W)
    rho_y <- sel_y$value
  }
  P_y <- solve_psd(crossprod(Theta) + rho_y * y_basis$W, t(Theta))
  
  blocks <- lapply(seq_len(P_number), function(p) {
    (X_centered[[p]] %*% t(P_x[[p]])) %*% x_basis$Psi
  })
  H <- do.call(cbind, c(blocks, if (Q_number > 0) list(w_centered) else NULL))
  D <- Y_centered %*% t(P_y)
  
  penalty_blocks <- c(rep(list(x_basis$W), P_number),
                      if (Q_number > 0) list(matrix(0, Q_number, Q_number))
                      else NULL)
  B <- block_diagonal(penalty_blocks)
  
  Omega <- crossprod(H)
  if (is.null(lambda)) {
    if (is.null(lambda_grid)) {
      lambda_grid <- default_penalty_grid(Omega, B)
    }
    scores <- numeric(length(lambda_grid))
    for (j in seq_along(lambda_grid)) {
      A_inverse_j <- inverse_psd(Omega + lambda_grid[j] * B)
      Gamma_j <- A_inverse_j %*% crossprod(H, D)
      effective_df <- sum(diag(A_inverse_j %*% Omega))
      denominator <- max(1 - effective_df / n, 1e-8)
      residual_coefficients <- D - H %*% Gamma_j
      integrated_error <- sum(
        (residual_coefficients %*% y_basis$Psi) * residual_coefficients
      ) / n
      scores[j] <- integrated_error / denominator^2
    }
    lambda <- lambda_grid[which.min(scores)]
    lambda_selection <- list(value = lambda, grid = lambda_grid, gcv = scores)
  } else {
    lambda_selection <- list(value = lambda, grid = lambda, gcv = NA_real_)
  }
  
  A <- symmetrize(Omega + lambda * B)
  A_inverse <- inverse_psd(A)
  Gamma_hat <- A_inverse %*% crossprod(H, D)        
  effective_df <- sum(diag(A_inverse %*% Omega))
  
  singular_values <- svd(H, nu = 0, nv = 0)$d
  rank_tolerance <- max(dim(H)) * max(singular_values) * .Machine$double.eps
  design_rank <- sum(singular_values > rank_tolerance)
  design_condition <- if (design_rank < L) {
    Inf
  } else {
    max(singular_values) / min(singular_values)
  }
  if (design_rank < L) {
    warning("The local mixed-predictor design is rank deficient (rank ",
            design_rank, " < ", L, "); ridge regularisation makes the fit ",
            "computable, but individual effects may not be identifiable")
  } else if (design_condition > 1e8) {
    warning("The local mixed-predictor design is ill-conditioned (condition ",
            "number ", format(design_condition, digits = 4), ")")
  }
  
  beta_on_grids <- lapply(seq_len(P_number), function(p) {
    rows <- ((p - 1) * Mx + 1):(p * Mx)
    Phi_x[[p]] %*% Gamma_hat[rows, , drop = FALSE] %*% t(Theta)   # Jx_p x Jy
  })
  EY <- matrix(0, n, length(y_grid))
  for (p in seq_len(P_number)) {
    J_p <- length(x_grids[[p]])
    paper_weights <- c(0.5, rep(1, max(J_p - 2, 0)), 0.5) / J_p
    EY <- EY + X_centered[[p]] %*% (paper_weights * beta_on_grids[[p]])
  }
  if (Q_number > 0) {
    rows <- (P_number * Mx + 1):L
    alpha_on_grid <- Gamma_hat[rows, , drop = FALSE] %*% t(Theta) 
    EY <- EY + w_centered %*% alpha_on_grid
  }
  ED <- EY %*% t(P_y)
  expected_gamma <- as.numeric(A_inverse %*% crossprod(H, ED))    
  
  dim_gamma <- L * My
  A_inv_H <- H %*% A_inverse
  if (variance_method == "HC") {
    if (n <= effective_df) {
      stop("HC variance requires n to exceed the effective degrees of freedom")
    }
    if (n < dim_gamma) {
      warning("n < L * My: the HC variance of vec(Gamma_hat) is rank ",
              "deficient; it is regularised, but consider variance_method ",
              "= 'model' or smaller bases")
    }
    correction <- n / (n - effective_df)
    R_conditional <- Y_centered - EY
    variance_gamma <- matrix(0, dim_gamma, dim_gamma)
    for (i in seq_len(n)) {
      u <- as.numeric(P_y %*% R_conditional[i, ])                 
      v <- A_inv_H[i, ]                                           
      variance_gamma <- variance_gamma +
        correction * kronecker(tcrossprod(u), tcrossprod(v))
    }
  } else {
    fitted_curves <- H %*% Gamma_hat %*% t(Theta)
    R_model <- Y_centered - fitted_curves
    response_covariance <- crossprod(R_model) / n                 
    measurement_covariance <- matrix(0, ncol(Y), ncol(Y))
    
    for (p in seq_len(P_number)) {
      J_p <- length(x_grids[[p]])
      weights_p <- c(0.5, rep(1, max(J_p - 2, 0)), 0.5) / J_p
      weighted_beta <- beta_on_grids[[p]] * weights_p
      measurement_covariance <- measurement_covariance +
        measurement_variances[p] * crossprod(weighted_beta)
    }
    
    conditional_response_covariance <-
      symmetrize(response_covariance + measurement_covariance)
    middle <- P_y %*% conditional_response_covariance %*% t(P_y)  # My x My
    outer_sum <- crossprod(A_inv_H)
    variance_gamma <- kronecker(middle, outer_sum)
  }
  variance_gamma <- regularize_psd(symmetrize(variance_gamma))
  
  structure(list(
    call = match.call(),
    gamma_hat = as.numeric(Gamma_hat),
    expected_gamma = expected_gamma,
    variance_gamma = variance_gamma,
    Gamma_matrix = Gamma_hat,
    y_mean_curve = y_mean_curve,
    x_mean_curves = x_mean_curves,
    w_means = w_means,
    y_grid = y_grid,
    x_grids = x_grids,
    P_x = P_x,
    Phi_x = Phi_x,
    P_y = P_y,
    Theta = Theta,
    x_basis = x_basis,
    y_basis = y_basis,
    rho_x = rho_x_used,
    rho_y = rho_y,
    lambda = lambda,
    lambda_selection = lambda_selection,
    variance_method = variance_method,
    measurement_variances = measurement_variances,
    effective_df = effective_df,
    design_rank = design_rank,
    design_condition = design_condition,
    predictor_names = predictor_names,
    covariate_names = covariate_names,
    basis_signature = basis_signature(
      x_basis, y_basis, predictor_names, covariate_names
    ),
    n = n,
    P_number = P_number,
    Q_number = Q_number,
    Mx = Mx,
    My = My,
    L = L
  ), class = "local_fofr_fit")
}
