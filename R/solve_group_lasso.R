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
