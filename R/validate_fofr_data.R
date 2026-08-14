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
