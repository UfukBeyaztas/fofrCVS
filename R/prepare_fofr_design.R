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
