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
    initial_basis <- bs(
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
