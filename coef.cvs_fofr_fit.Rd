\name{coef.cvs_fofr_fit}
\alias{coef.cvs_fofr_fit}
\title{Extract Basis Coefficients from a CVS-FoFR Fit}
\description{
Extracts the tensor-product basis coefficient matrix associated with the
Target-only, control-variates (CVS), or penalized control-variates (pCVS)
estimator from a fitted function-on-function transfer model.
}
\usage{
\method{coef}{cvs_fofr_fit}(
  object,
  method = c("local", "CVS", "pCVS"),
  ...
)
}
\arguments{
  \item{object}{
An object of class \code{"cvs_fofr_fit"}, as returned by
\code{\link{fit_cvs_fofr}} or \code{\link{combine_cvs_fofr}}.
  }
  \item{method}{
A character string specifying the estimator to extract. \code{"local"}
returns the Target-only estimator, \code{"CVS"} returns the feasible
control-variates estimator, and \code{"pCVS"} returns the selected penalized
control-variates estimator. Partial matching is supported.
  }
  \item{...}{
Additional arguments for the \code{coef} generic. They are currently unused.
  }
}
\details{
Let \eqn{P} denote the number of functional predictors, \eqn{Q} the number of
scalar covariates, and \eqn{M_x} and \eqn{M_y} the predictor- and
response-basis dimensions. The fitted coefficient object is represented by a
matrix \eqn{\widehat{\boldsymbol{\Gamma}}} of dimension
\eqn{L \times M_y}, where \eqn{L=P M_x+Q}. Its first \eqn{P M_x} rows form
the functional blocks
\eqn{\widehat{\boldsymbol{\Gamma}}_p \in
\mathbb{R}^{M_x \times M_y}}, and its final \eqn{Q} rows form the scalar block
\eqn{\widehat{\boldsymbol{\Gamma}}_w \in
\mathbb{R}^{Q \times M_y}}.

For common predictor and response bases
\eqn{\boldsymbol{\phi}(s)} and \eqn{\boldsymbol{\vartheta}(t)}, respectively,
the fitted coefficient functions are reconstructed as
\deqn{
\widehat\beta_p(s,t)
=\boldsymbol{\phi}(s)^\mathsf{T}
 \widehat{\boldsymbol{\Gamma}}_p
 \boldsymbol{\vartheta}(t),
}
and
\deqn{
\widehat{\boldsymbol{\alpha}}(t)
=\widehat{\boldsymbol{\Gamma}}_w\boldsymbol{\vartheta}(t).
}
Thus, the matrix returned by this method contains coefficients in the common
finite-dimensional sieve; it does not contain evaluations of the coefficient
functions on a grid.

For \code{method = "local"}, the matrix is the penalized Target-only
estimator. For \code{method = "CVS"}, it is obtained by reshaping the feasible
control-variates coefficient vector
\deqn{
\widehat{\boldsymbol{g}}_{C}
=\widehat{\boldsymbol{g}}_0
-\widehat{\boldsymbol{U}}_{\star}
 \left(\widehat{\boldsymbol{\Delta}}-\widehat{\boldsymbol{d}}\right),
}
where the adjustment combines centered target--source coefficient differences
using the estimated study-level covariance matrices. For
\code{method = "pCVS"}, the returned matrix corresponds to the selected point
on the pCVS path. The latter estimates the source-specific discrepancies by a
group lasso after transformation by the square root of the coefficient Gram
matrix. Consequently, the group penalty measures integrated coefficient
distance and is invariant to rescaling of the basis.

The stored estimator is vectorized and reshaped columnwise, consistently with
the convention \eqn{\boldsymbol{g}=\mathrm{vec}(\boldsymbol{\Gamma})}
used in the transfer construction. Extraction with \code{method = "pCVS"}
requires that a pCVS solution has been selected, either at a user-supplied
penalty or by a target validation sample; otherwise the method returns an
error.

The numerical entries of \eqn{\widehat{\boldsymbol{\Gamma}}} depend on the
chosen bases. Interpretation and comparison of fitted effects should therefore
usually be based on \code{\link{coefficient_surface}} or
\code{\link{scalar_effect}}, which reconstruct the corresponding functions.
The row names retain the fitted predictor and covariate names so that the
block structure can be identified without consulting the original data.
}
\value{
A numeric matrix with \eqn{P M_x+Q} rows and \eqn{M_y} columns. For each
functional predictor, \eqn{M_x} consecutive rows contain its predictor-basis
coefficients; the final \eqn{Q} rows correspond to the scalar covariates. The
columns index the response-basis functions. The matrix has an attribute named
\code{"method"} recording the extracted estimator.
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)

target <- simulate_fofr_data(40, s_grid = s_grid, t_grid = t_grid)
source <- simulate_fofr_data(80, s_grid = s_grid, t_grid = t_grid)
validation <- simulate_fofr_data(40, s_grid = s_grid, t_grid = t_grid)

fit <- fit_cvs_fofr(
  datasets = list(target, source),
  x_basis_type = "bspline",
  y_basis_type = "bspline",
  Mx = 5,
  My = 5,
  variance_method = "model",
  target_validation = validation
)

gamma_target <- coef(fit, method = "local")
gamma_cvs <- coef(fit, method = "CVS")
gamma_pcvs <- coef(fit, method = "pCVS")

dim(gamma_pcvs)
head(gamma_pcvs)
attr(gamma_pcvs, "method")
}
}
