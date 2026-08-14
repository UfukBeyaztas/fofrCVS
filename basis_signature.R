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