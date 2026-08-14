\name{print.cvs_fofr_fit}
\alias{print.cvs_fofr_fit}
\title{Print a CVS-FoFR Fit}
\description{
Prints a compact summary of a fitted function-on-function control-variates
transfer model, including model dimensions, target-site regularization and
design diagnostics, and the status of pCVS penalty selection.
}
\usage{
\method{print}{cvs_fofr_fit}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"cvs_fofr_fit"}, as returned by
\code{\link{fit_cvs_fofr}} or \code{\link{combine_cvs_fofr}}.
  }
  \item{...}{
Additional arguments for the \code{print} generic. They are currently unused.
  }
}
\details{
The method reports the target sample size; the numbers of functional
predictors, scalar covariates, and source studies; and the predictor- and
response-basis dimensions \eqn{M_x} and \eqn{M_y}. If there are \eqn{P}
functional predictors and \eqn{Q} scalar covariates, the mixed-predictor
design has
\deqn{
L=P M_x+Q
}
columns, while the vectorized coefficient object has dimension \eqn{L M_y}.

The displayed target value of \code{lambda} is the roughness penalty used in
the local target fit. It is either the user-supplied value or the value chosen
by generalized cross-validation in \code{\link{fit_local_fofr}}. The target
design rank is the numerical rank of the centered mixed-predictor design
matrix, reported relative to its \eqn{L} columns. The condition number is the
ratio of the largest to the smallest retained singular value when the design
has full column rank; it is reported as infinite for a rank-deficient design.
These quantities diagnose numerical conditioning and identifiability of the
local coefficient blocks. They are not measures of predictive performance.

If no pCVS penalty has been selected, the method prints
\sQuote{pCVS: no zeta selected}. This occurs when neither a user-specified
\code{zeta} nor a target validation sample was supplied to the combining
procedure. In that case, the Target-only and CVS estimators and the fitted
pCVS path remain available, but there is no selected pCVS coefficient vector
for extraction or prediction.

When a pCVS estimator has been selected, the method reports its penalty
\eqn{\zeta}, the selection rule, and the relative
Karush--Kuhn--Tucker (KKT) residual of the selected path solution. The pCVS
groups are regularized in the Gram geometry of the coefficient functions. If
\eqn{\boldsymbol{\Psi}_x} and \eqn{\boldsymbol{\Psi}_y} are the predictor and
response Gram matrices, this geometry is induced by
\deqn{
\boldsymbol{G}
=\boldsymbol{\Psi}_y\otimes
\mathrm{diag}(\boldsymbol{\Psi}_x,\ldots,
\boldsymbol{\Psi}_x,\boldsymbol{I}_Q).
}
Consequently, the group norm corresponds to the integrated \eqn{L^2} norm of
the collection of coefficient surfaces and scalar-effect functions rather
than an unscaled Euclidean norm of basis coefficients. The relative KKT
residual is a dimensionless numerical optimality diagnostic; smaller values
indicate closer satisfaction of the first-order conditions. It should not be
interpreted as an estimation error, a test statistic, or a measure of model
fit.

The output is intentionally selective. Coefficient estimates, validation
losses, the complete pCVS path, convergence records, and source-specific
discrepancies are stored in \code{x} but are not printed by this method.
}
\value{
The object \code{x}, returned invisibly. The method is called for its printed
side effect.
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)

target <- simulate_fofr_data(35, s_grid = s_grid, t_grid = t_grid)
source <- simulate_fofr_data(60, s_grid = s_grid, t_grid = t_grid)
validation <- simulate_fofr_data(35, s_grid = s_grid, t_grid = t_grid)

fit <- fit_cvs_fofr(
  datasets = list(target, source),
  x_basis_type = "bspline",
  y_basis_type = "bspline",
  Mx = 5,
  My = 5,
  variance_method = "model",
  target_validation = validation
)

print(fit)

returned_object <- print(fit)
identical(returned_object, fit)
}
}
