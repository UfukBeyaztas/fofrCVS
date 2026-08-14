\name{scalar_effect}
\alias{scalar_effect}
\title{Evaluate a Scalar Effect from a CVS-FoFR Fit}
\description{
Reconstructs the response-varying coefficient function associated with one
scalar covariate from the Target-only, control-variates (CVS), or penalized
control-variates (pCVS) estimator in a fitted function-on-function transfer
model.
}
\usage{
scalar_effect(
  object,
  covariate = 1,
  method = c("local", "CVS", "pCVS"),
  t_grid = NULL
)
}
\arguments{
  \item{object}{
An object of class \code{"cvs_fofr_fit"}, as returned by
\code{\link{fit_cvs_fofr}} or \code{\link{combine_cvs_fofr}}.
  }
  \item{covariate}{
The scalar covariate whose coefficient function is to be evaluated. Supply
either a single integer between 1 and the number of scalar covariates or the
exact covariate name recorded in \code{object}.
  }
  \item{method}{
A character string specifying the estimator. \code{"local"} uses the
Target-only estimator, \code{"CVS"} uses the feasible control-variates
estimator, and \code{"pCVS"} uses the selected penalized control-variates
estimator. Partial matching is supported.
  }
  \item{t_grid}{
An optional numeric vector of response-domain locations at which the scalar
effect is evaluated. If \code{NULL}, 201 equally spaced points over the
boundary of the fitted response basis are used.
  }
}
\details{
For a function-on-function regression model with \eqn{P} functional
predictors and \eqn{Q} scalar covariates, the conditional mean has the form
\deqn{
\mathbb{E}\{Y_i(t)\mid X_{i1},\ldots,X_{iP},\boldsymbol{w}_i\}
=\mu_0(t)
+\sum_{p=1}^{P}\int_{\mathcal{S}}
 X_{ip}(s)\beta_p(s,t)\,ds
+\sum_{q=1}^{Q}w_{iq}\alpha_q(t).
}
The function \eqn{\alpha_q(t)} is therefore the change in the conditional
mean response at location \eqn{t} associated with a one-unit increase in
scalar covariate \eqn{q}, holding the other predictors fixed. The fitted model
is estimated after within-study centering, but centering changes the functional
intercept rather than the coefficient function. Hence the same slope
interpretation applies relative to the target covariate mean.

Let \eqn{M_x} and \eqn{M_y} denote the predictor- and response-basis
dimensions and put \eqn{L=P M_x+Q}. The fitted coefficient matrix
\eqn{\widehat{\boldsymbol{\Gamma}}} has \eqn{L} rows. Its first \eqn{P M_x}
rows contain the coefficient blocks of the functional predictors, and its
final \eqn{Q} rows contain the scalar-effect coefficients. If
\eqn{\widehat{\boldsymbol{\gamma}}_{wq}^{\mathsf{T}}} is the row associated
with covariate \eqn{q} and \eqn{\boldsymbol{\vartheta}(t)} is the response
basis, the reconstructed effect is
\deqn{
\widehat{\alpha}_q(t)
=\widehat{\boldsymbol{\gamma}}_{wq}^{\mathsf{T}}
 \boldsymbol{\vartheta}(t).
}
The function selects the requested coefficient vector, reshapes it columnwise
into \eqn{\widehat{\boldsymbol{\Gamma}}}, extracts row
\eqn{P M_x+q}, and evaluates this expression on \code{t_grid}.

For \code{method = "local"}, the result is obtained from the penalized
Target-only fit. For \code{method = "CVS"}, it uses the feasible
control-variates adjustment based on the centered target--source coefficient
differences. For \code{method = "pCVS"}, it uses the selected point on the
group-penalized transfer path. The pCVS discrepancy penalty acts jointly on
all coefficient surfaces and scalar-effect functions contributed by a source;
it does not select scalar effects separately. Within the coefficient Gram
metric, the scalar part uses the response Gram matrix, so its squared norm is
the integrated squared magnitude of \eqn{\alpha_q(t)}.

The supplied \code{t_grid} is an evaluation grid only: changing it does not
refit the model. For interpretation on the fitted response domain, its values
should lie within the boundary of the response basis. If the model contains no
scalar covariates, the function returns an error. Extraction with
\code{method = "pCVS"} requires that a pCVS penalty has been selected by a
target validation sample or supplied through \code{zeta}; otherwise an error
is returned.
}
\value{
A data frame with one row per response-domain evaluation location and the
following columns:
  \item{t_grid}{The response-domain evaluation locations.}
  \item{alpha}{The reconstructed scalar coefficient function.}
  \item{method}{The estimator used for reconstruction.}
  \item{covariate}{The numerical index of the selected scalar covariate.}
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

continuous_effect <- scalar_effect(
  fit,
  covariate = "continuous",
  method = "pCVS",
  t_grid = t_grid
)

head(continuous_effect)

plot(
  t_grid,
  target$truth$alpha$continuous,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = expression(alpha[1](t))
)
lines(
  continuous_effect$t_grid,
  continuous_effect$alpha,
  col = "red",
  lwd = 2
)
legend(
  "topright",
  legend = c("Truth", "pCVS"),
  col = c("black", "red"),
  lwd = 2,
  bty = "n"
)
}
}
