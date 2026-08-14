\name{fit_local_fofr}
\alias{fit_local_fofr}
\title{Fit a Local Function-on-Function Regression Model}
\description{
Fits the study-specific component of the control-variates transfer procedure.
The function reconstructs the observed curves in common predictor and response
bases, estimates the mixed-predictor function-on-function regression model by
penalized least squares, and constructs the coefficient mean and covariance
summaries required for CVS and pCVS transfer.
}
\usage{
fit_local_fofr(
  dataset,
  x_basis,
  y_basis,
  rho_x = NULL,
  rho_y = NULL,
  lambda = NULL,
  lambda_grid = NULL,
  variance_method = c("HC", "model")
)
}
\arguments{
  \item{dataset}{
A study-level list containing \code{Y}, \code{y_grid}, \code{X_list}, and
\code{x_grids}, and optionally a scalar-covariate matrix \code{w}; see
Details.
  }
  \item{x_basis}{
A predictor-domain basis object created by \code{\link{fofr_basis}}.
  }
  \item{y_basis}{
A response-domain basis object created by \code{\link{fofr_basis}}.
  }
  \item{rho_x}{
Predictor-curve smoothing penalty. Supply \code{NULL}, one nonnegative number
used for every functional predictor, or a nonnegative vector with one value
per functional predictor. If \code{NULL}, a penalty is selected separately
for each predictor by generalized cross-validation.
  }
  \item{rho_y}{
Response-curve smoothing penalty. Supply a single nonnegative number or
\code{NULL}. If \code{NULL}, it is selected by generalized
cross-validation.
  }
  \item{lambda}{
A single nonnegative roughness penalty for the function-on-function
coefficient surfaces. If \code{NULL}, the value is selected by generalized
cross-validation over \code{lambda_grid}.
  }
  \item{lambda_grid}{
An optional nonempty numeric vector of nonnegative candidates for
\code{lambda}. It is used only when \code{lambda = NULL}. If both are
\code{NULL}, a data-scaled logarithmic grid is constructed internally.
  }
  \item{variance_method}{
Method used to estimate the conditional covariance of the local coefficient
vector. \code{"HC"} gives a heteroskedasticity-consistent estimator based on
subject-specific residual outer products. \code{"model"} uses a common
residual covariance augmented by a correction for predictor-reconstruction
noise. Partial matching is supported.
  }
}
\details{
The response \code{Y} must be an \eqn{n \times J_y} matrix whose rows are
subjects and whose columns correspond to the strictly increasing vector
\code{y_grid}. The component \code{X_list} must contain at least one
functional-predictor matrix. Its \eqn{p}th element has \eqn{n} rows and one
column for every location in the strictly increasing vector
\code{x_grids[[p]]}. The optional component \code{w} is an \eqn{n \times Q}
matrix of scalar covariates. All observations and grid values must be finite,
and at least five subjects are required. Names attached to \code{X_list} and
the columns of \code{w} are retained in the fitted object; otherwise default
names are assigned.

Let \eqn{\boldsymbol{\phi}(s)} and
\eqn{\boldsymbol{\vartheta}(t)} denote the predictor and response bases, with
dimensions \eqn{M_x} and \eqn{M_y}. Their Gram matrices are
\eqn{\boldsymbol{\Psi}_x} and \eqn{\boldsymbol{\Psi}_y}, and their integrated
second-derivative penalty matrices are \eqn{\boldsymbol{W}_x} and
\eqn{\boldsymbol{W}_y}. Every predictor grid must lie inside the boundary of
\code{x_basis}, and the response grid must lie inside that of \code{y_basis}.

The response, functional predictors, and scalar covariates are centered within
the study. For predictor \eqn{p}, let \eqn{\boldsymbol{\Phi}_p} be the
evaluation matrix of \eqn{\boldsymbol{\phi}} on \code{x_grids[[p]]}. If
\eqn{\boldsymbol{x}_{ip}} is the centered observed trajectory, its smoothed
basis coefficients are
\deqn{
\widehat{\boldsymbol{c}}_{ip}
=\boldsymbol{P}_{xp}\boldsymbol{x}_{ip},
\qquad
\boldsymbol{P}_{xp}
=(\boldsymbol{\Phi}_p^\mathsf{T}\boldsymbol{\Phi}_p
+\rho_{xp}\boldsymbol{W}_x)^{-1}
\boldsymbol{\Phi}_p^\mathsf{T}.
}
Similarly, if \eqn{\boldsymbol{\Theta}} evaluates
\eqn{\boldsymbol{\vartheta}} on \code{y_grid}, the smoothed response
coefficients are
\deqn{
\widehat{\boldsymbol{\xi}}_i
=\boldsymbol{P}_y\boldsymbol{y}_i,
\qquad
\boldsymbol{P}_y
=(\boldsymbol{\Theta}^\mathsf{T}\boldsymbol{\Theta}
+\rho_y\boldsymbol{W}_y)^{-1}
\boldsymbol{\Theta}^\mathsf{T}.
}
When a smoothing penalty is not supplied, it minimizes the usual
generalized-cross-validation score based on the corresponding linear
smoother.

Put \eqn{L=P M_x+Q}. The integrated design vector for subject \eqn{i} is
\deqn{
\boldsymbol{h}_i
=\left(\widehat{\boldsymbol{c}}_{i1}^\mathsf{T}
\boldsymbol{\Psi}_x,\ldots,
\widehat{\boldsymbol{c}}_{iP}^\mathsf{T}
\boldsymbol{\Psi}_x,\boldsymbol{w}_i^\mathsf{T}\right)^\mathsf{T},
}
where \eqn{\boldsymbol{w}_i} is centered. Let \eqn{\boldsymbol{H}} have rows
\eqn{\boldsymbol{h}_i^\mathsf{T}} and let \eqn{\boldsymbol{D}} have rows
\eqn{\widehat{\boldsymbol{\xi}}_i^\mathsf{T}}. The coefficient matrix
\eqn{\boldsymbol{\Gamma}\in\mathbb{R}^{L\times M_y}} is partitioned into
\eqn{P} functional blocks followed by a \eqn{Q\times M_y} scalar block. The
local criterion is
\deqn{
\mathrm{tr}\{(\boldsymbol{D}-\boldsymbol{H}\boldsymbol{\Gamma})
\boldsymbol{\Psi}_y
(\boldsymbol{D}-\boldsymbol{H}\boldsymbol{\Gamma})^\mathsf{T}\}
+\lambda\,
\mathrm{tr}(\boldsymbol{\Gamma}^\mathsf{T}\boldsymbol{B}
\boldsymbol{\Gamma}\boldsymbol{\Psi}_y),
}
where
\eqn{\boldsymbol{B}=\mathrm{diag}(\boldsymbol{W}_x,\ldots,
\boldsymbol{W}_x,\boldsymbol{0}_{Q\times Q})}. Thus, roughness in the
predictor direction is penalized for every coefficient surface, while the
scalar-effect functions are unpenalized in this step. With
\eqn{\boldsymbol{\Omega}=\boldsymbol{H}^\mathsf{T}\boldsymbol{H}} and
\eqn{\boldsymbol{A}=\boldsymbol{\Omega}+\lambda\boldsymbol{B}}, the fitted
matrix is
\deqn{
\widehat{\boldsymbol{\Gamma}}
=\boldsymbol{A}^{-1}\boldsymbol{H}^\mathsf{T}\boldsymbol{D}.
}
If \code{lambda} is not supplied, generalized cross-validation uses the
integrated residual coefficient loss and effective degrees of freedom
\eqn{\nu=\mathrm{tr}(\boldsymbol{A}^{-1}\boldsymbol{\Omega})}.

In exact arithmetic, uniqueness of the penalized solution requires
\deqn{
\mathrm{Null}(\boldsymbol{H})\cap
\mathrm{Null}(\lambda\boldsymbol{B})=\{\boldsymbol{0}\}.
}
The implementation reports the numerical rank and condition number of
\eqn{\boldsymbol{H}} and warns about deficient or severely ill-conditioned
designs. Spectrally floored inverses keep such fits computable, but individual
coefficient effects may remain weakly identified.

The transfer summary requires an estimate of the conditional mean and
covariance of
\eqn{\widehat{\boldsymbol{g}}=
\mathrm{vec}(\widehat{\boldsymbol{\Gamma}})}. The component
\code{expected_gamma} is a study-computable plug-in mean: the fitted
conditional response curves are reconstructed on the observation grid,
smoothed through \eqn{\boldsymbol{P}_y}, and mapped back through the local
normal equations. It should not in general be interpreted as the exact
conditional expectation of \code{gamma_hat}; exactness requires additional
working-model and tuning conditions.

For \code{variance_method = "HC"}, the covariance estimate has the form
\deqn{
\widehat{\boldsymbol{\Sigma}}_{HC}
=\frac{n}{n-\nu}\sum_{i=1}^{n}
(\widehat{\boldsymbol{r}}_i\widehat{\boldsymbol{r}}_i^\mathsf{T})
\otimes
(\boldsymbol{a}_i\boldsymbol{a}_i^\mathsf{T}),
\qquad
\boldsymbol{a}_i=\boldsymbol{A}^{-1}\boldsymbol{h}_i,
}
where \eqn{\widehat{\boldsymbol{r}}_i} is the smoothed residual coefficient
vector. This option requires \eqn{n>\nu}; before spectral regularization it is
rank deficient when \eqn{n<L M_y}. For \code{variance_method = "model"}, a
common response covariance, augmented by the estimated predictor-reconstruction
contribution, is combined with
\eqn{\boldsymbol{A}^{-1}\boldsymbol{H}^\mathsf{T}
\boldsymbol{H}\boldsymbol{A}^{-1}} through a Kronecker product. When tuning
parameters are selected using the same responses employed in fitting, these
covariance formulas are working plug-in estimates rather than exact
conditional identities.
}
\value{
An object of class \code{"local_fofr_fit"}. It is a list containing:
  \item{call}{The matched function call.}
  \item{gamma_hat}{
The columnwise vectorization of the fitted coefficient matrix, of length
\eqn{(P M_x+Q)M_y}.
  }
  \item{expected_gamma}{
The plug-in conditional-mean vector used in the feasible transfer summary.
  }
  \item{variance_gamma}{
The spectrally regularized estimated covariance matrix of \code{gamma_hat}.
  }
  \item{Gamma_matrix}{
The fitted \eqn{(P M_x+Q)\times M_y} coefficient matrix.
  }
  \item{y_mean_curve}{The study-level mean response on \code{y_grid}.}
  \item{x_mean_curves}{
A list containing the study-level mean curve of every functional predictor.
  }
  \item{w_means}{
The scalar-covariate means, or a numeric vector of length zero when the model
has no scalar covariates.
  }
  \item{y_grid}{The response observation grid.}
  \item{x_grids}{The functional-predictor observation grids.}
  \item{P_x}{
A list of predictor smoothing matrices, one for each functional predictor.
  }
  \item{Phi_x}{
A list of predictor-basis evaluation matrices on the observed predictor grids.
  }
  \item{P_y}{The response smoothing matrix.}
  \item{Theta}{The response-basis evaluation matrix on \code{y_grid}.}
  \item{x_basis}{The predictor basis supplied in \code{x_basis}.}
  \item{y_basis}{The response basis supplied in \code{y_basis}.}
  \item{rho_x}{The predictor smoothing penalties actually used.}
  \item{rho_y}{The response smoothing penalty actually used.}
  \item{lambda}{The coefficient-surface roughness penalty actually used.}
  \item{lambda_selection}{
A list with components \code{value}, \code{grid}, and \code{gcv}. When
\code{lambda} was supplied, \code{grid} contains that value and \code{gcv}
is \code{NA_real_}.
  }
  \item{variance_method}{The covariance estimation method.}
  \item{measurement_variances}{
The mean squared predictor-reconstruction residual for each functional
predictor.
  }
  \item{effective_df}{
The effective degrees of freedom
\eqn{\mathrm{tr}(\boldsymbol{A}^{-1}\boldsymbol{\Omega})}.
  }
  \item{design_rank}{The numerical rank of the mixed-predictor design matrix.}
  \item{design_condition}{
The condition number of the design matrix, or \code{Inf} when it is rank
deficient.
  }
  \item{predictor_names}{The ordered functional-predictor names.}
  \item{covariate_names}{The ordered scalar-covariate names.}
  \item{basis_signature}{
A compact record of the basis types, dimensions, boundaries, knots, Gram and
penalty matrices, and ordered predictor and covariate names. It is used to
verify compatibility across studies.
  }
  \item{n}{The number of subjects in the study.}
  \item{P_number}{The number of functional predictors.}
  \item{Q_number}{The number of scalar covariates.}
  \item{Mx}{The predictor-basis dimension.}
  \item{My}{The response-basis dimension.}
  \item{L}{The number of rows in the coefficient matrix, \eqn{P M_x+Q}.}
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)
study_data <- simulate_fofr_data(50, s_grid = s_grid, t_grid = t_grid)

x_basis <- fofr_basis("bspline", M = 5)
y_basis <- fofr_basis("bspline", M = 5)

local_fit <- fit_local_fofr(
  dataset = study_data,
  x_basis = x_basis,
  y_basis = y_basis,
  variance_method = "model"
)

dim(local_fit$Gamma_matrix)
local_fit$lambda
local_fit$effective_df
c(rank = local_fit$design_rank, columns = local_fit$L)

# This is the object transmitted by a source site.
site_summary <- fofr_site_summary(local_fit)
names(site_summary)
}
}
