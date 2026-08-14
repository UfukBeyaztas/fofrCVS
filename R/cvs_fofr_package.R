






solve_psd <- function(A, b, relative_floor = 1e-10) {
  inverse_psd(A, relative_floor = relative_floor) %*% b
}

trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2 || any(diff(grid) <= 0)) {
    stop("grid must contain at least two strictly increasing values")
  }
  
  increments <- diff(grid)
  weights <- numeric(length(grid))
  weights[1] <- increments[1] / 2
  weights[length(grid)] <- increments[length(increments)] / 2
  
  if (length(grid) > 2) {
    weights[2:(length(grid) - 1)] <-
      (increments[1:(length(increments) - 1)] + increments[2:length(increments)]) / 2
  }
  
  weights
}

block_diagonal <- function(matrices) {
  dimensions <- vapply(matrices, nrow, numeric(1))
  answer <- matrix(0, sum(dimensions), sum(dimensions))
  starts <- cumsum(c(1, dimensions[-length(dimensions)]))
  
  for (j in seq_along(matrices)) {
    indices <- starts[j]:(starts[j] + dimensions[j] - 1)
    answer[indices, indices] <- matrices[[j]]
  }
  
  answer
}

fofr_basis <- function(type = c("fourier", "bspline"),
                       M = 15,
                       degree = 3,
                       boundary = c(0, 1),
                       integration_grid = seq(boundary[1], boundary[2],
                                              length.out = 1001)) {
  type <- match.arg(type)

  if (length(M) != 1 || !is.finite(M) || M != as.integer(M)) {
    stop("M must be a single integer")
  }
  if (length(degree) != 1 || !is.finite(degree) ||
      degree != as.integer(degree) || degree < 2) {
    stop("degree must be a single integer of at least 2")
  }
  if (length(boundary) != 2 || any(!is.finite(boundary)) ||
      boundary[1] >= boundary[2]) {
    stop("boundary must contain two increasing finite values")
  }
  integration_grid <- as.numeric(integration_grid)
  if (length(integration_grid) < 3 || any(!is.finite(integration_grid)) ||
      any(diff(integration_grid) <= 0) ||
      min(integration_grid) < boundary[1] ||
      max(integration_grid) > boundary[2]) {
    stop("integration_grid must contain increasing finite values within boundary")
  }
  
  if (M < 3) {
    stop("M must be at least 3")
  }
  if (type == "bspline" && M < degree + 1) {
    stop("For a B-spline basis, M must be at least degree + 1")
  }
  
  basis <- list(
    type = type,
    M = M,
    degree = degree,
    boundary = boundary,
    integration_grid = integration_grid
  )
  
  if (type == "bspline") {
    initial_basis <- splines::bs(
      integration_grid,
      df = M,
      degree = degree,
      intercept = TRUE,
      Boundary.knots = boundary
    )
    basis$knots <- attr(initial_basis, "knots")
    basis$full_knots <- c(
      rep(boundary[1], degree + 1),
      basis$knots,
      rep(boundary[2], degree + 1)
    )
  }
  
  Phi <- evaluate_basis(basis, integration_grid)
  Phi_second <- evaluate_basis_second_derivative(basis, integration_grid)
  integration_weights <- trapezoid_weights(integration_grid)
  
  basis$Psi <- symmetrize(crossprod(Phi, Phi * integration_weights))
  basis$W <- symmetrize(crossprod(Phi_second, Phi_second * integration_weights))
  structure(basis, class = "fofr_basis")
}

evaluate_basis <- function(basis, grid) {
  grid <- as.numeric(grid)
  
  if (basis$type == "bspline") {
    return(splines::splineDesign(
      knots = basis$full_knots,
      x = grid,
      ord = basis$degree + 1,
      derivs = rep(0, length(grid)),
      outer.ok = TRUE
    ))
  }
  
  domain_length <- basis$boundary[2] - basis$boundary[1]
  u <- (grid - basis$boundary[1]) / domain_length
  
  Phi <- matrix(0, nrow = length(grid), ncol = basis$M)
  Phi[, 1] <- 1 / sqrt(domain_length)
  
  for (m in 2:basis$M) {
    if (m %% 2 == 0) {
      frequency <- m / 2
      Phi[, m] <- sqrt(2 / domain_length) * cos(2 * pi * frequency * u)
    } else {
      frequency <- (m - 1) / 2
      Phi[, m] <- sqrt(2 / domain_length) * sin(2 * pi * frequency * u)
    }
  }
  
  Phi
}

evaluate_basis_second_derivative <- function(basis, grid) {
  grid <- as.numeric(grid)
  
  if (basis$type == "bspline") {
    return(splines::splineDesign(
      knots = basis$full_knots,
      x = grid,
      ord = basis$degree + 1,
      derivs = rep(2, length(grid)),
      outer.ok = TRUE
    ))
  }
  
  domain_length <- basis$boundary[2] - basis$boundary[1]
  u <- (grid - basis$boundary[1]) / domain_length
  
  Phi_second <- matrix(0, nrow = length(grid), ncol = basis$M)
  
  for (m in 2:basis$M) {
    if (m %% 2 == 0) {
      frequency <- m / 2
      Phi_second[, m] <- -(2 * pi * frequency / domain_length)^2 *
        sqrt(2 / domain_length) * cos(2 * pi * frequency * u)
    } else {
      frequency <- (m - 1) / 2
      Phi_second[, m] <- -(2 * pi * frequency / domain_length)^2 *
        sqrt(2 / domain_length) * sin(2 * pi * frequency * u)
    }
  }
  
  Phi_second
}

default_penalty_grid <- function(data_matrix, penalty_matrix,
                                 lower_power = -6, upper_power = 2,
                                 length.out = 31) {
  data_scale <- sum(diag(data_matrix)) / nrow(data_matrix)
  penalty_scale <- sum(diag(penalty_matrix)) / nrow(penalty_matrix)
  
  if (!is.finite(penalty_scale) || penalty_scale <= .Machine$double.eps) {
    penalty_scale <- 1
  }
  
  reference <- max(data_scale / penalty_scale, .Machine$double.eps)
  reference * 10^seq(lower_power, upper_power, length.out = length.out)
}

select_curve_smoothing <- function(Z, Phi, W, rho_grid = NULL) {
  if (is.null(rho_grid)) {
    rho_grid <- default_penalty_grid(crossprod(Phi), W)
  }
  
  n <- nrow(Z)
  J <- ncol(Z)
  scores <- numeric(length(rho_grid))
  
  for (j in seq_along(rho_grid)) {
    P <- solve_psd(crossprod(Phi) + rho_grid[j] * W, t(Phi))
    smoother <- Phi %*% P
    residuals <- Z - Z %*% t(smoother)
    effective_df <- sum(diag(smoother))
    denominator <- max(1 - effective_df / J, 1e-8)
    scores[j] <- mean(residuals^2) / denominator^2
  }
  
  list(
    value = rho_grid[which.min(scores)],
    grid = rho_grid,
    gcv = scores
  )
}

# ---------------------------------------------------------------------------
# pCVS optimization utilities
# ---------------------------------------------------------------------------

solve_group_lasso <- function(delta_hat,
                                       precision_delta,
                                       zeta,
                                       group_size,
                                       number_groups,
                                       initial = NULL,
                                       tolerance = 1e-8,
                                       kkt_tolerance = 1e-6,
                                       max_iterations = 20000,
                                       lipschitz_constant = NULL) {
  if (length(zeta) != 1 || !is.finite(zeta) || zeta < 0) {
    stop("zeta must be a single nonnegative finite number")
  }
  if (length(delta_hat) != group_size * number_groups) {
    stop("length(delta_hat) must equal group_size * number_groups")
  }
  if (length(tolerance) != 1 || !is.finite(tolerance) || tolerance <= 0 ||
      length(kkt_tolerance) != 1 || !is.finite(kkt_tolerance) ||
      kkt_tolerance <= 0 || length(max_iterations) != 1 ||
      !is.finite(max_iterations) || max_iterations < 1) {
    stop("tolerances and max_iterations must be positive")
  }
  
  if (zeta == 0) {
    return(list(
      delta = as.numeric(delta_hat),
      iterations = 0L,
      converged = TRUE,
      kkt_residual = 0,
      relative_kkt_residual = 0
    ))
  }
  
  precision_delta <- (as.matrix(precision_delta) +
                        t(as.matrix(precision_delta))) / 2
  if (is.null(lipschitz_constant)) {
    largest_eigenvalue <- max(eigen(
      precision_delta, symmetric = TRUE, only.values = TRUE
    )$values)
    lipschitz_constant <- max(2 * largest_eigenvalue, 1e-12)
  } else if (length(lipschitz_constant) != 1 ||
             !is.finite(lipschitz_constant) || lipschitz_constant <= 0) {
    stop("lipschitz_constant must be a positive finite number")
  }
  step_size <- 1 / lipschitz_constant
  
  group_threshold <- function(value, threshold) {
    answer <- value
    for (k in seq_len(number_groups)) {
      indices <- ((k - 1) * group_size + 1):(k * group_size)
      group <- value[indices]
      group_norm <- sqrt(sum(group^2))
      answer[indices] <- if (group_norm <= threshold) {
        0
      } else {
        (1 - threshold / group_norm) * group
      }
    }
    answer
  }
  
  delta <- if (is.null(initial)) {
    rep(0, length(delta_hat))
  } else {
    as.numeric(initial)
  }
  if (length(delta) != length(delta_hat) || any(!is.finite(delta))) {
    stop("initial must be finite and have the same length as delta_hat")
  }
  
  accelerated <- delta
  acceleration <- 1
  converged <- FALSE
  kkt_residual <- Inf
  relative_kkt_residual <- Inf
  kkt_scale <- max(
    1,
    zeta,
    max(abs(2 * as.numeric(precision_delta %*% delta_hat)))
  )
  
  for (iteration in seq_len(max_iterations)) {
    gradient <- 2 * as.numeric(
      precision_delta %*% (accelerated - delta_hat)
    )
    proposed <- group_threshold(
      accelerated - step_size * gradient,
      step_size * zeta
    )
    
    next_acceleration <- (1 + sqrt(1 + 4 * acceleration^2)) / 2
    next_accelerated <- proposed +
      ((acceleration - 1) / next_acceleration) * (proposed - delta)
    if (sum((accelerated - proposed) * (proposed - delta)) > 0) {
      next_acceleration <- 1
      next_accelerated <- proposed
    }
    
    relative_change <- sqrt(sum((proposed - delta)^2)) /
      max(1, sqrt(sum(delta^2)))
    proposed_gradient <- 2 * as.numeric(
      precision_delta %*% (proposed - delta_hat)
    )
    kkt_residual <- 0
    
    for (k in seq_len(number_groups)) {
      indices <- ((k - 1) * group_size + 1):(k * group_size)
      group <- proposed[indices]
      group_norm <- sqrt(sum(group^2))
      if (group_norm > 1e-12) {
        group_residual <- proposed_gradient[indices] +
          zeta * group / group_norm
        kkt_residual <- max(
          kkt_residual,
          sqrt(sum(group_residual^2))
        )
      } else {
        kkt_residual <- max(
          kkt_residual,
          max(0, sqrt(sum(proposed_gradient[indices]^2)) - zeta)
        )
      }
    }
    relative_kkt_residual <- kkt_residual / kkt_scale
    
    delta <- proposed
    accelerated <- next_accelerated
    acceleration <- next_acceleration
    if (relative_change < tolerance &&
        relative_kkt_residual < kkt_tolerance) {
      converged <- TRUE
      break
    }
  }
  
  list(
    delta = delta,
    iterations = iteration,
    converged = converged,
    kkt_residual = kkt_residual,
    relative_kkt_residual = relative_kkt_residual
  )
}

penalty_path <- function(delta_hat, precision_delta,
                           group_size, number_groups,
                           min_ratio = 1e-3, length.out = 31) {
  gradient_at_zero <- 2 * as.numeric(precision_delta %*% delta_hat)
  group_norms <- numeric(number_groups)
  
  for (k in seq_len(number_groups)) {
    indices <- ((k - 1) * group_size + 1):(k * group_size)
    group_norms[k] <- sqrt(sum(gradient_at_zero[indices]^2))
  }
  
  zeta_max <- max(group_norms)
  
  if (!is.finite(zeta_max) || zeta_max <= .Machine$double.eps) {
    return(0)
  }
  
  exp(seq(log(zeta_max), log(zeta_max * min_ratio), length.out = length.out))
}

# ---------------------------------------------------------------------------
# Function-on-function regression and transfer-learning functions
# ---------------------------------------------------------------------------

validate_fofr_data <- function(dataset) {
  required <- c("Y", "y_grid", "X_list", "x_grids")
  if (!all(required %in% names(dataset))) {
    stop("Each dataset must contain Y, y_grid, X_list, and x_grids ",
         "(and optionally a scalar covariate matrix w)")
  }
  
  Y <- as.matrix(dataset$Y)
  y_grid <- as.numeric(dataset$y_grid)
  X_list <- lapply(dataset$X_list, as.matrix)
  x_grids <- lapply(dataset$x_grids, as.numeric)
  n <- nrow(Y)
  
  if (ncol(Y) != length(y_grid)) {
    stop("ncol(Y) must equal length(y_grid)")
  }
  if (any(!is.finite(y_grid))) {
    stop("y_grid must contain only finite values")
  }
  if (length(X_list) != length(x_grids)) {
    stop("X_list and x_grids must have the same length")
  }
  if (length(X_list) < 1) {
    stop("At least one functional predictor is required")
  }
  for (p in seq_along(X_list)) {
    if (nrow(X_list[[p]]) != n) {
      stop("Each functional predictor must have one row per subject")
    }
    if (ncol(X_list[[p]]) != length(x_grids[[p]])) {
      stop("Functional predictor ", p, ": ncol must equal its grid length")
    }
    if (any(diff(x_grids[[p]]) <= 0)) {
      stop("Functional predictor grids must be strictly increasing")
    }
    if (any(!is.finite(x_grids[[p]]))) {
      stop("Functional predictor grids must contain only finite values")
    }
    if (any(!is.finite(X_list[[p]]))) {
      stop("Functional predictors must contain only finite values")
    }
  }
  if (any(diff(y_grid) <= 0)) {
    stop("y_grid must be strictly increasing")
  }
  if (any(!is.finite(Y))) {
    stop("Y must contain only finite values")
  }
  
  w <- dataset$w
  if (!is.null(w)) {
    w <- as.matrix(w)
    if (nrow(w) != n) {
      stop("w must have one row per subject")
    }
    if (any(!is.finite(w))) {
      stop("w must contain only finite values")
    }
  }
  if (n < 5) {
    stop("At least five subjects are required in each dataset")
  }
  
  predictor_names <- names(X_list)
  if (is.null(predictor_names) || any(predictor_names == "")) {
    predictor_names <- paste0("X", seq_along(X_list))
  }
  names(X_list) <- predictor_names
  names(x_grids) <- predictor_names
  
  covariate_names <- if (is.null(w)) character(0) else colnames(w)
  if (!is.null(w) && (is.null(covariate_names) || any(covariate_names == ""))) {
    covariate_names <- paste0("w", seq_len(ncol(w)))
    colnames(w) <- covariate_names
  }
  
  list(
    Y = Y,
    y_grid = y_grid,
    X_list = X_list,
    x_grids = x_grids,
    w = w,
    predictor_names = predictor_names,
    covariate_names = covariate_names
  )
}

matrix_root <- function(A, inverse = FALSE,
                               relative_floor = 1e-10,
                               absolute_floor = 1e-12) {
  A <- symmetrize(as.matrix(A))
  eig <- eigen(A, symmetric = TRUE)
  scale <- max(abs(eig$values))
  floor_value <- max(absolute_floor, relative_floor * scale)
  values <- pmax(eig$values, floor_value)
  weights <- if (inverse) 1 / sqrt(values) else sqrt(values)
  symmetrize(eig$vectors %*% (t(eig$vectors) * weights))
}

basis_signature <- function(x_basis, y_basis, predictor_names,
                                 covariate_names) {
  compact <- function(basis) {
    list(
      type = basis$type,
      M = basis$M,
      degree = basis$degree,
      boundary = basis$boundary,
      knots = basis$knots,
      Psi = basis$Psi,
      W = basis$W
    )
  }
  
  list(
    x_basis = compact(x_basis),
    y_basis = compact(y_basis),
    predictor_names = predictor_names,
    covariate_names = covariate_names
  )
}

# Local fit: two-step ridge mirroring Section 2 of the paper.
# Step 1 smooths each predictor curve onto the x-basis (penalty rho_x) and the
# response curves onto the y-basis (penalty rho_y). Step 2 solves
#   Gamma_hat = argmin tr{(Y P_y' - H Gamma) Psi_y
#                         (Y P_y' - H Gamma)'}
#                     + lambda tr(Gamma' B Gamma Psi_y)
# in closed form, where H stacks the predictor blocks [Psi_x b_i1', ..., w_i']
# and B penalises roughness of each beta_p(., t) in the s direction. Because
# Psi_y is positive definite, this has the same normal equation as the
# unweighted Frobenius criterion, while GCV below uses integrated loss.
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
  
  # ---- centering ----
  y_mean_curve <- colMeans(Y)
  Y_centered <- sweep(Y, 2, y_mean_curve, "-")
  x_mean_curves <- lapply(X_list, colMeans)
  X_centered <- lapply(seq_len(P_number), function(p) {
    sweep(X_list[[p]], 2, x_mean_curves[[p]], "-")
  })
  w_means <- if (Q_number > 0) colMeans(w) else numeric(0)
  w_centered <- if (Q_number > 0) sweep(w, 2, w_means, "-") else NULL
  
  # ---- step 1: smoothing ----
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
  
  # ---- design matrix H (n x L) and smoothed response D (n x My) ----
  blocks <- lapply(seq_len(P_number), function(p) {
    (X_centered[[p]] %*% t(P_x[[p]])) %*% x_basis$Psi
  })
  H <- do.call(cbind, c(blocks, if (Q_number > 0) list(w_centered) else NULL))
  D <- Y_centered %*% t(P_y)
  
  # roughness penalty in the s direction, block-diagonal over predictors,
  # zero block for scalar covariates
  penalty_blocks <- c(rep(list(x_basis$W), P_number),
                      if (Q_number > 0) list(matrix(0, Q_number, Q_number))
                      else NULL)
  B <- block_diagonal(penalty_blocks)
  
  # ---- step 2: ridge with GCV over lambda ----
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
  Gamma_hat <- A_inverse %*% crossprod(H, D)        # L x My
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
  
  # ---- conditional mean of Y given Z (plug-in, Appendix A.2 analogue) ----
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
    alpha_on_grid <- Gamma_hat[rows, , drop = FALSE] %*% t(Theta) # Q x Jy
    EY <- EY + w_centered %*% alpha_on_grid
  }
  ED <- EY %*% t(P_y)
  expected_gamma <- as.numeric(A_inverse %*% crossprod(H, ED))    # vec, col-major
  
  # ---- conditional variance of vec(Gamma_hat) given Z ----
  # vec(Gamma_hat) = sum_i (P_y y_i) (x) (A^{-1} h_i), so with rows of Y
  # independent given Z,
  #   Var(vec Gamma_hat | Z) = sum_i (P_y Sigma_i P_y') (x)
  #                                   (A^{-1} h_i h_i' A^{-1}).
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
      u <- as.numeric(P_y %*% R_conditional[i, ])                  # My
      v <- A_inv_H[i, ]                                            # L
      variance_gamma <- variance_gamma +
        correction * kronecker(tcrossprod(u), tcrossprod(v))
    }
  } else {
    fitted_curves <- H %*% Gamma_hat %*% t(Theta)
    R_model <- Y_centered - fitted_curves
    response_covariance <- crossprod(R_model) / n                  # Jy x Jy
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

fofr_site_summary <- function(local_fit) {
  if (!inherits(local_fit, "local_fofr_fit")) {
    stop("local_fit must be created by fit_local_fofr")
  }
  structure(list(
    gamma_hat = local_fit$gamma_hat,
    expected_gamma = local_fit$expected_gamma,
    variance_gamma = local_fit$variance_gamma,
    Mx = local_fit$Mx,
    My = local_fit$My,
    P_number = local_fit$P_number,
    Q_number = local_fit$Q_number,
    basis_signature = local_fit$basis_signature
  ), class = "cvs_fofr_summary")
}

as_fofr_site_summary <- function(object) {
  if (inherits(object, "local_fofr_fit")) {
    return(fofr_site_summary(object))
  }
  if (!inherits(object, "cvs_fofr_summary")) {
    stop("Sources must be local_fofr_fit or cvs_fofr_summary objects")
  }
  object
}

coefficient_metric <- function(local_fit) {
  row_blocks <- c(
    rep(list(local_fit$x_basis$Psi), local_fit$P_number),
    if (local_fit$Q_number > 0) list(diag(local_fit$Q_number)) else NULL
  )
  row_metric <- block_diagonal(row_blocks)
  regularize_psd(symmetrize(kronecker(local_fit$y_basis$Psi, row_metric)))
}

prepare_fofr_design <- function(target_fit, X_list, w = NULL,
                                x_grids = NULL) {
  P_number <- target_fit$P_number
  Q_number <- target_fit$Q_number
  if (!is.list(X_list)) {
    X_list <- list(X_list)
  }
  if (length(X_list) != P_number) {
    stop("newdata must contain ", P_number, " functional predictors")
  }
  input_predictor_names <- names(X_list)
  if (!is.null(input_predictor_names) &&
      all(nzchar(input_predictor_names))) {
    if (!setequal(input_predictor_names, target_fit$predictor_names)) {
      stop("Named functional predictors must match the target predictor names")
    }
    X_list <- X_list[target_fit$predictor_names]
  }
  if (is.null(x_grids)) {
    x_grids <- target_fit$x_grids
  }
  if (!is.list(x_grids)) {
    x_grids <- list(x_grids)
  }
  if (length(x_grids) != P_number) {
    stop("x_grids must contain one grid per functional predictor")
  }
  input_grid_names <- names(x_grids)
  if (!is.null(input_grid_names) && all(nzchar(input_grid_names))) {
    if (!setequal(input_grid_names, target_fit$predictor_names)) {
      stop("Named x_grids must match the target predictor names")
    }
    x_grids <- x_grids[target_fit$predictor_names]
  }
  x_grids <- lapply(x_grids, as.numeric)
  
  X_list <- lapply(seq_len(P_number), function(p) {
    Xp <- X_list[[p]]
    if (is.null(dim(Xp))) {
      Xp <- matrix(Xp, nrow = 1)
    }
    Xp <- as.matrix(Xp)
    if (ncol(Xp) != length(x_grids[[p]])) {
      stop("Functional predictor ", p,
           ": ncol must equal the length of its new-data grid")
    }
    if (any(!is.finite(Xp)) || any(!is.finite(x_grids[[p]])) ||
        any(diff(x_grids[[p]]) <= 0)) {
      stop("New functional data and grids must be finite, with increasing grids")
    }
    Xp
  })
  n_new <- nrow(X_list[[1]])
  if (any(vapply(X_list, nrow, integer(1)) != n_new)) {
    stop("All new functional predictors must have the same number of rows")
  }
  
  blocks <- lapply(seq_len(P_number), function(p) {
    new_grid <- x_grids[[p]]
    training_grid <- target_fit$x_grids[[p]]
    grid_tolerance <- 100 * .Machine$double.eps *
      max(1, max(abs(training_grid)))
    if (min(new_grid) < min(training_grid) - grid_tolerance ||
        max(new_grid) > max(training_grid) + grid_tolerance) {
      stop("Functional predictor ", p,
           " uses a grid outside the target training-grid range")
    }
    
    same_grid <- length(new_grid) == length(training_grid) &&
      max(abs(new_grid - training_grid)) <= grid_tolerance
    if (same_grid) {
      mean_curve <- target_fit$x_mean_curves[[p]]
      P_new <- target_fit$P_x[[p]]
    } else {
      mean_curve <- approx(training_grid, target_fit$x_mean_curves[[p]],
                           xout = new_grid, rule = 2)$y
      Phi_new <- evaluate_basis(target_fit$x_basis, new_grid)
      P_new <- solve_psd(
        crossprod(Phi_new) + target_fit$rho_x[p] * target_fit$x_basis$W,
        t(Phi_new)
      )
    }
    Xc <- sweep(X_list[[p]], 2, mean_curve, "-")
    (Xc %*% t(P_new)) %*% target_fit$x_basis$Psi
  })
  
  if (Q_number > 0) {
    if (is.null(w)) {
      stop("newdata must include the scalar covariates w")
    }
    if (is.null(dim(w))) {
      w <- matrix(w, nrow = n_new)
    }
    w <- as.matrix(w)
    if (ncol(w) != Q_number || nrow(w) != n_new) {
      stop("w must be ", n_new, " x ", Q_number)
    }
    if (!is.null(colnames(w)) && all(nzchar(colnames(w)))) {
      if (!setequal(colnames(w), target_fit$covariate_names)) {
        stop("Named scalar covariates must match the target covariate names")
      }
      w <- w[, target_fit$covariate_names, drop = FALSE]
    }
    if (any(!is.finite(w))) {
      stop("w must contain only finite values")
    }
    blocks <- c(blocks, list(sweep(w, 2, target_fit$w_means, "-")))
  }
  do.call(cbind, blocks)
}

prepare_response_grid <- function(target_fit, y_grid = NULL) {
  if (is.null(y_grid)) {
    y_grid <- target_fit$y_grid
  }
  y_grid <- as.numeric(y_grid)
  if (length(y_grid) < 2 || any(!is.finite(y_grid)) ||
      any(diff(y_grid) <= 0)) {
    stop("The response grid must be finite and strictly increasing")
  }
  training_grid <- target_fit$y_grid
  grid_tolerance <- 100 * .Machine$double.eps *
    max(1, max(abs(training_grid)))
  if (min(y_grid) < min(training_grid) - grid_tolerance ||
      max(y_grid) > max(training_grid) + grid_tolerance) {
    stop("The response grid is outside the target training-grid range")
  }
  same_grid <- length(y_grid) == length(training_grid) &&
    max(abs(y_grid - training_grid)) <= grid_tolerance
  mean_curve <- if (same_grid) {
    target_fit$y_mean_curve
  } else {
    approx(training_grid, target_fit$y_mean_curve,
           xout = y_grid, rule = 2)$y
  }
  list(
    grid = y_grid,
    basis_matrix = evaluate_basis(target_fit$y_basis, y_grid),
    mean_curve = mean_curve
  )
}

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
  
  # pCVS must be invariant to arbitrary basis scaling.  In theta coordinates,
  # Euclidean group norms equal the integrated L2 norms of all coefficient
  # surfaces and scalar-effect functions:
  #   ||theta||^2 = vec(Gamma)' (Psi_y (x) R_x) vec(Gamma).
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
  
  # Include the exact target-only endpoint.
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

# Generate data from the mixed-predictor function-on-function model used in
# the package examples. Setting effect = "shifted" changes the regression
# effects and therefore produces a nontransferable source. A covariance shift
# is obtained by retaining effect = "common" and changing covariance_scale.
simulate_fofr_data <- function(n,
                               s_grid = seq(0, 1, length.out = 101),
                               t_grid = seq(0, 1, length.out = 101),
                               covariance_scale = 4,
                               effect = c("common", "shifted"),
                               measurement_sd = 0.05,
                               error_sd = 0.10) {
  effect <- match.arg(effect)
  if (length(n) != 1 || !is.finite(n) || n != as.integer(n) || n < 5) {
    stop("n must be a single integer of at least 5")
  }
  s_grid <- as.numeric(s_grid)
  t_grid <- as.numeric(t_grid)
  if (length(s_grid) < 2 || any(!is.finite(s_grid)) ||
      any(diff(s_grid) <= 0)) {
    stop("s_grid must contain at least two increasing finite values")
  }
  if (length(t_grid) < 2 || any(!is.finite(t_grid)) ||
      any(diff(t_grid) <= 0)) {
    stop("t_grid must contain at least two increasing finite values")
  }
  if (length(covariance_scale) != 1 || !is.finite(covariance_scale) ||
      covariance_scale <= 0) {
    stop("covariance_scale must be a positive finite number")
  }
  if (length(measurement_sd) != 1 || !is.finite(measurement_sd) ||
      measurement_sd < 0 || length(error_sd) != 1 ||
      !is.finite(error_sd) || error_sd < 0) {
    stop("measurement_sd and error_sd must be nonnegative finite numbers")
  }

  number_s <- length(s_grid)
  number_t <- length(t_grid)
  covariance_matrix <- covariance_scale *
    exp(-10 * abs(outer(s_grid, s_grid, "-")))
  covariance_factor <- chol(
    covariance_matrix + diag(1e-9, number_s)
  )

  latent_activity <- matrix(rnorm(n * number_s), nrow = n) %*%
    covariance_factor
  latent_exposure <- matrix(rnorm(n * number_s), nrow = n) %*%
    covariance_factor
  w <- cbind(
    continuous = rnorm(n),
    binary = rbinom(n, size = 1, prob = 0.5)
  )

  beta_activity <- outer(
    s_grid, t_grid,
    function(s, t) sin(2 * pi * s) * cos(pi * t)
  )
  beta_exposure <- outer(
    s_grid, t_grid,
    function(s, t) 2 * (s - 0.5) * (t - 0.5)
  )
  alpha_continuous <- 1 + 0.5 * sin(2 * pi * t_grid)
  alpha_binary <- -0.8 * t_grid

  if (effect == "shifted") {
    beta_activity <- -0.35 * beta_activity
    beta_exposure <- beta_exposure + outer(
      s_grid, t_grid,
      function(s, t) 1.25 * cos(pi * s) * sin(2 * pi * t)
    )
    alpha_continuous <- -0.50 * alpha_continuous
    alpha_binary <- alpha_binary + 0.9 * sin(pi * t_grid)
  }

  integration_weights <- trapezoid_weights(s_grid)
  conditional_mean <-
    latent_activity %*% (integration_weights * beta_activity) +
    latent_exposure %*% (integration_weights * beta_exposure) +
    w[, "continuous"] %o% alpha_continuous +
    w[, "binary"] %o% alpha_binary
  Y <- conditional_mean +
    matrix(rnorm(n * number_t, sd = error_sd), nrow = n)

  observed_activity <- latent_activity +
    matrix(rnorm(n * number_s, sd = measurement_sd), nrow = n)
  observed_exposure <- latent_exposure +
    matrix(rnorm(n * number_s, sd = measurement_sd), nrow = n)

  list(
    Y = Y,
    y_grid = t_grid,
    X_list = list(
      activity = observed_activity,
      exposure = observed_exposure
    ),
    x_grids = list(activity = s_grid, exposure = s_grid),
    w = w,
    conditional_mean = conditional_mean,
    latent_X_list = list(
      activity = latent_activity,
      exposure = latent_exposure
    ),
    truth = list(
      beta = list(
        activity = beta_activity,
        exposure = beta_exposure
      ),
      alpha = list(
        continuous = alpha_continuous,
        binary = alpha_binary
      )
    ),
    effect = effect,
    covariance_scale = covariance_scale
  )
}

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
