\name{fit_cvs_fofr}
\alias{fit_cvs_fofr}
\title{Fit CVS and pCVS Function-on-Function Transfer Models}
\description{
Fits the complete control-variates transfer procedure for function-on-function
linear regression with multiple functional predictors and optional scalar
covariates. The first data set is treated as the target study and all remaining
data sets as sources. The returned object contains the Target-only, feasible
control-variates (CVS), and penalized control-variates (pCVS) estimators.
}
\usage{
fit_cvs_fofr(
  datasets,
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
  max_iterations = 20000
)
}
\arguments{
  \item{datasets}{
A list containing the target data set first and at least one source data set
thereafter. Each data set must contain \code{Y}, \code{y_grid},
\code{X_list}, and \code{x_grids}, and may contain a scalar-covariate matrix
\code{w}; see Details.
  }
  \item{x_basis_type}{
The basis used over the functional-predictor domain: \code{"fourier"} or
\code{"bspline"}. Partial matching is supported.
  }
  \item{y_basis_type}{
The basis used over the response domain: \code{"fourier"} or
\code{"bspline"}. Partial matching is supported.
  }
  \item{Mx}{
The number of predictor-basis functions. It must be at least 3 and, for a
B-spline basis, at least \code{degree + 1}.
  }
  \item{My}{
The number of response-basis functions. It must be at least 3 and, for a
B-spline basis, at least \code{degree + 1}.
  }
  \item{degree}{
The polynomial degree used for both B-spline bases. It must be an integer of
at least 2. This argument does not affect a Fourier basis.
  }
  \item{x_boundary}{
A length-two increasing numeric vector defining the common
functional-predictor domain. Every predictor grid in every study must lie
within this interval.
  }
  \item{y_boundary}{
A length-two increasing numeric vector defining the common response domain.
Every response grid in every study must lie within this interval.
  }
  \item{rho_x}{
Predictor-curve smoothing penalty. Supply \code{NULL}, one nonnegative number
used for all functional predictors, or a nonnegative vector with one value per
functional predictor. If \code{NULL}, a penalty is selected separately for
each predictor in each study by generalized cross-validation.
  }
  \item{rho_y}{
Response-curve smoothing penalty. Supply a single nonnegative number or
\code{NULL}. If \code{NULL}, it is selected separately in each study by
generalized cross-validation.
  }
  \item{lambda}{
The nonnegative roughness penalty for the function-on-function coefficient
surfaces. If supplied, the same value is used in every local study. If
\code{NULL}, each study selects its own value by generalized
cross-validation over \code{lambda_grid}.
  }
  \item{lambda_grid}{
An optional nonempty numeric vector of nonnegative candidate values for
\code{lambda}. It is used only when \code{lambda = NULL}. If both arguments
are \code{NULL}, a data-scaled logarithmic grid is constructed separately for
each study.
  }
  \item{variance_method}{
Method used to estimate the conditional covariance of each local coefficient
vector. \code{"HC"} gives a heteroskedasticity-consistent residual
outer-product estimator with an effective-degrees-of-freedom correction.
\code{"model"} uses a common residual covariance augmented by a correction
for functional-predictor reconstruction noise. Partial matching is supported.
  }
  \item{zeta}{
An optional nonnegative finite pCVS penalty. It is added to the candidate path
and is selected when \code{target_validation = NULL}. If validation data are
provided, it is treated as an additional candidate and validation selection
takes precedence.
  }
  \item{zeta_grid}{
An optional numeric vector of nonnegative finite pCVS penalties. If
\code{NULL}, a decreasing data-adaptive path is constructed. The exact
Target-only endpoint \code{zeta = 0} is always included.
  }
  \item{target_validation}{
An optional independent target validation data set with the same component
structure as an element of \code{datasets}. It is used only to select the
pCVS penalty by integrated squared prediction error and is not included in any
local fit.
  }
  \item{tolerance}{
A positive finite convergence tolerance for the accelerated proximal-gradient
algorithm used to compute the pCVS path.
  }
  \item{max_iterations}{
A positive number giving the maximum number of proximal-gradient iterations
at each pCVS penalty.
  }
}
\details{
The function accommodates the centered mixed-predictor model
\deqn{
Y_i^{(k)}(t)-\mu_Y^{(k)}(t)
=\sum_{p=1}^{P}\int_{\mathcal{S}}
\{X_{ip}^{(k)}(s)-\mu_{X_p}^{(k)}(s)\}
\beta_p^{(k)}(s,t)\,ds
+\{W_i^{(k)}-\mu_W^{(k)}\}^\mathsf{T}
\boldsymbol{\alpha}^{(k)}(t)+\varepsilon_i^{(k)}(t),
}
where study 0 is the target and studies \eqn{1,\ldots,K} are sources. The
coefficient surfaces and scalar-effect functions may differ between the target
and the sources; equality of the regression operators is not imposed.

Each element of \code{datasets} is a study-level list. The response
\code{Y} is an \eqn{n_k \times J_y^{(k)}} matrix whose rows are subjects and
whose columns correspond to the strictly increasing vector \code{y_grid}.
The component \code{X_list} is a list of \eqn{P} matrices; its \eqn{p}th
matrix has \eqn{n_k} rows and one column for each location in
\code{x_grids[[p]]}. The optional \code{w} is an \eqn{n_k \times Q} matrix.
All entries must be finite and each study must contain at least five subjects.
The numbers and ordered names of the functional predictors and scalar
covariates must agree across studies. Observation grids may differ by study
and predictor, provided that they remain inside the common basis boundaries.

The same predictor and response bases are constructed once and used in every
study. Within each study, the recorded predictor and response curves are first
smoothed onto these bases. If \eqn{\boldsymbol{\Phi}_p} and
\eqn{\boldsymbol{\Theta}} are the corresponding evaluation matrices, the
coefficient smoothers have the form
\deqn{
\boldsymbol{P}_{xp}
=(\boldsymbol{\Phi}_p^\mathsf{T}\boldsymbol{\Phi}_p
+\rho_{xp}\boldsymbol{W}_x)^{-1}
\boldsymbol{\Phi}_p^\mathsf{T},
\qquad
\boldsymbol{P}_y
=(\boldsymbol{\Theta}^\mathsf{T}\boldsymbol{\Theta}
+\rho_y\boldsymbol{W}_y)^{-1}
\boldsymbol{\Theta}^\mathsf{T}.
}
The resulting local coefficient matrix minimizes
\deqn{
\mathrm{tr}\{(\boldsymbol{D}_k-\boldsymbol{H}_k\boldsymbol{\Gamma})
\boldsymbol{\Psi}_y
(\boldsymbol{D}_k-\boldsymbol{H}_k\boldsymbol{\Gamma})^\mathsf{T}\}
+\lambda_k\,
\mathrm{tr}(\boldsymbol{\Gamma}^\mathsf{T}\boldsymbol{B}
\boldsymbol{\Gamma}\boldsymbol{\Psi}_y),
}
where \eqn{\boldsymbol{B}} applies the second-derivative penalty to the
functional blocks and leaves the scalar-effect block unpenalized. Writing
\eqn{\boldsymbol{A}_k=\boldsymbol{H}_k^\mathsf{T}\boldsymbol{H}_k+
\lambda_k\boldsymbol{B}}, the fitted matrix is
\deqn{
\widehat{\boldsymbol{\Gamma}}_k
=\boldsymbol{A}_k^{-1}\boldsymbol{H}_k^\mathsf{T}\boldsymbol{D}_k.
}
The implementation uses spectrally stabilized inverses when a fitted matrix is
numerically singular or nearly singular and records local rank and condition
diagnostics in \code{local_fits}.

After fitting all studies, the function extracts each local coefficient
estimate, plug-in conditional mean, covariance estimate, and basis signature,
then calls \code{\link{combine_cvs_fofr}}. The Target-only estimator is the
coefficient estimate from \code{datasets[[1]]}. CVS subtracts an optimally
weighted centered target--source difference. pCVS estimates source-specific
mean discrepancies by a group lasso after transformation by the coefficient
Gram matrix. In that metric, the penalty is the integrated Hilbert--Schmidt
distance between coefficient objects and does not depend on arbitrary scaling
of the bases.

The pCVS path always contains \code{zeta = 0}, at which the estimator is
exactly Target-only. When \code{target_validation} is supplied, only converged
path solutions satisfying the Karush--Kuhn--Tucker tolerance are eligible, and
the value minimizing integrated validation curve MSE is selected. When no
validation sample is supplied, \code{zeta} must be specified to select a pCVS
estimator. If both \code{target_validation} and \code{zeta} are \code{NULL},
the path is returned but \code{gamma_pcvs} is \code{NULL}.
}
\value{
An object of class \code{"cvs_fofr_fit"}. It contains all components
documented for \code{\link{combine_cvs_fofr}}, including the Target-only,
CVS, and selected pCVS coefficient vectors; the pCVS coefficient and
discrepancy paths; the coefficient Gram metric and control-variate weights;
validation losses; and convergence, iteration, and KKT diagnostics. In
addition, this interface returns:
  \item{local_fits}{
A list of \code{"local_fofr_fit"} objects in the same order as
\code{datasets}. The first element is the target fit and the remaining
elements are the source fits.
  }
  \item{call}{The matched function call.}
The principal coefficient components are \code{gamma_local},
\code{gamma_cvs}, and \code{gamma_pcvs}. The last is \code{NULL} when no
pCVS penalty was selected.
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)

target <- simulate_fofr_data(35, s_grid = s_grid, t_grid = t_grid)
source_1 <- simulate_fofr_data(60, s_grid = s_grid, t_grid = t_grid)
source_2 <- simulate_fofr_data(
  60,
  s_grid = s_grid,
  t_grid = t_grid,
  effect = "shifted"
)
validation <- simulate_fofr_data(35, s_grid = s_grid, t_grid = t_grid)
new_target <- simulate_fofr_data(10, s_grid = s_grid, t_grid = t_grid)

fit <- fit_cvs_fofr(
  datasets = list(target, source_1, source_2),
  x_basis_type = "bspline",
  y_basis_type = "bspline",
  Mx = 5,
  My = 5,
  variance_method = "model",
  target_validation = validation
)

fit
fit$selected_zeta
coef(fit, method = "pCVS")[1:4, 1:3]

predicted_curves <- predict(fit, new_target, method = "pCVS")
dim(predicted_curves)

activity_surface <- coefficient_surface(
  fit,
  predictor = "activity",
  method = "pCVS"
)
dim(activity_surface$beta)
}
}
