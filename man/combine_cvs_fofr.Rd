\name{combine_cvs_fofr}
\alias{combine_cvs_fofr}
\title{Combine Local Fits by CVS and pCVS Transfer}
\description{
Combines a target-site function-on-function fit with one or more source-site
summaries. The function computes the Target-only and feasible
control-variates (CVS) estimators, constructs the penalized control-variates
(pCVS) solution path, and optionally selects a pCVS estimator using an
independent target validation sample or a user-specified penalty.
}
\usage{
combine_cvs_fofr(
  target_fit,
  source_summaries,
  zeta = NULL,
  zeta_grid = NULL,
  target_validation = NULL,
  tolerance = 1e-7,
  max_iterations = 20000
)
}
\arguments{
  \item{target_fit}{
An object of class \code{"local_fofr_fit"} obtained by applying
\code{\link{fit_local_fofr}} to the target study.
  }
  \item{source_summaries}{
A nonempty list of source-site summaries produced by
\code{\link{fofr_site_summary}}. Objects of class \code{"local_fofr_fit"}
are also accepted and are converted internally. All studies must use
identical predictor and response bases and the same ordered functional
predictors and scalar covariates.
  }
  \item{zeta}{
An optional nonnegative finite number specifying a pCVS penalty. The value is
added to the candidate path. If \code{target_validation = NULL}, the pCVS
estimate at this value is selected. If validation data are supplied,
\code{zeta} is treated as an additional candidate and validation selection
takes precedence.
  }
  \item{zeta_grid}{
An optional numeric vector of nonnegative finite pCVS penalties. If
\code{NULL}, a decreasing data-adaptive path is constructed from the
source-discrepancy gradients. The exact Target-only endpoint
\code{zeta = 0} is always included, duplicate values are removed, and the
resulting path is sorted in decreasing order.
  }
  \item{target_validation}{
An optional independent target validation data set. It must be a list with
components \code{Y}, \code{y_grid}, \code{X_list}, and \code{x_grids}, and
with \code{w} when scalar covariates are present. The pCVS penalty is selected
by integrated squared prediction error on these response curves.
  }
  \item{tolerance}{
A positive finite convergence tolerance for the accelerated proximal-gradient
algorithm used along the pCVS path.
  }
  \item{max_iterations}{
A positive number giving the maximum number of proximal-gradient iterations
at each value of \code{zeta}.
  }
}
\details{
Each study is represented by a coefficient vector
\eqn{\widehat{\boldsymbol{g}}_k}, a plug-in conditional mean
\eqn{\widehat{\boldsymbol{\mu}}_k}, an estimated covariance matrix
\eqn{\widehat{\boldsymbol{\Sigma}}_k}, and a basis signature. Let study 0
denote the target and let studies \eqn{1,\ldots,K} denote the sources. The
stacked target--source differences are
\deqn{
\widehat{\boldsymbol{\Delta}}
=
\left\lbrace
(\widehat{\boldsymbol{g}}_0-\widehat{\boldsymbol{g}}_1)^\mathsf{T},
\ldots,
(\widehat{\boldsymbol{g}}_0-\widehat{\boldsymbol{g}}_K)^\mathsf{T}
\right\rbrace^\mathsf{T},
}
with estimated mean
\eqn{\widehat{\boldsymbol{d}}} formed analogously from the exported plug-in
means. Writing
\deqn{
\widehat{\boldsymbol{\Sigma}}_\star
=\left(\sum_{k=0}^{K}
\widehat{\boldsymbol{\Sigma}}_k^{-1}\right)^{-1}
}
and
\deqn{
\widehat{\boldsymbol{U}}_\star
=\widehat{\boldsymbol{\Sigma}}_\star
\left[\widehat{\boldsymbol{\Sigma}}_1^{-1},\ldots,
\widehat{\boldsymbol{\Sigma}}_K^{-1}\right],
}
the feasible control-variates estimator is
\deqn{
\widehat{\boldsymbol{g}}_C
=\widehat{\boldsymbol{g}}_0-
\widehat{\boldsymbol{U}}_\star
(\widehat{\boldsymbol{\Delta}}-\widehat{\boldsymbol{d}}).
}
This construction combines precision-weighted, centered local fluctuations;
it does not pool subject-level observations or require the source regression
operators to equal the target operator.

The pCVS estimator regularizes the source discrepancies in the functional
geometry of the coefficient object. Let
\eqn{\boldsymbol{G}=\boldsymbol{\Psi}_y\otimes
\mathrm{diag}(\boldsymbol{\Psi}_x,\ldots,
\boldsymbol{\Psi}_x,\boldsymbol{I}_Q)} be the coefficient Gram matrix and
define metric coordinates
\eqn{\widehat{\boldsymbol{\theta}}_k=
\boldsymbol{G}^{1/2}\widehat{\boldsymbol{g}}_k}. For a stacked candidate
discrepancy
\eqn{\boldsymbol{u}=(\boldsymbol{u}_1^\mathsf{T},\ldots,
\boldsymbol{u}_K^\mathsf{T})^\mathsf{T}}, the path solves
\deqn{
\widehat{\boldsymbol{d}}_\zeta
=\mathop{\mathrm{arg\,min}}_{\boldsymbol{u}}
\left\{(\widehat{\boldsymbol{\Delta}}_\theta-\boldsymbol{u})^\mathsf{T}
\boldsymbol{Q}_{\mathrm{alg}}
(\widehat{\boldsymbol{\Delta}}_\theta-\boldsymbol{u})
+\zeta\sum_{k=1}^{K}\|\boldsymbol{u}_k\|_2\right\}.
}
The corresponding coefficient estimator is
\deqn{
\widehat{\boldsymbol{\theta}}_{PC}(\zeta)
=\widehat{\boldsymbol{\theta}}_0-
\widehat{\boldsymbol{U}}_{\theta\star}
(\widehat{\boldsymbol{\Delta}}_\theta-
\widehat{\boldsymbol{d}}_\zeta).
}
Because multiplication by \eqn{\boldsymbol{G}^{1/2}} is an isometry, each
group norm is the integrated Hilbert--Schmidt distance between the target and
one source coefficient object, rather than a basis-dependent Euclidean norm.

At \code{zeta = 0}, the discrepancy estimate equals the observed discrepancy
and pCVS reduces exactly to the Target-only estimator. At sufficiently large
penalties, the discrepancy blocks are zero and the sources enter through the
zero-discrepancy control-variate combination. Between these endpoints, a
source with a discrepancy block near zero transfers strongly, whereas a block
retained near its observed target--source difference largely withdraws that
source. The entire path is computed with warm starts. A solution is eligible
for validation selection only if both the convergence and Karush--Kuhn--Tucker
criteria are satisfied.

When \code{target_validation} is supplied, the selected penalty minimizes
\deqn{
\frac{1}{n_v}\sum_{i=1}^{n_v}
\int_{\mathcal{T}}
\{Y_i(t)-\widehat{Y}_{i,\zeta}(t)\}^2\,dt
}
over the converged path solutions, using trapezoidal integration on the
validation response grid. If neither \code{target_validation} nor
\code{zeta} is supplied, the CVS estimator and the full pCVS path are returned
but no pCVS coefficient vector is selected; consequently,
\code{gamma_pcvs} and \code{selected_zeta} are \code{NULL}.

Before combination, the function verifies dimensions, basis definitions,
Gram and penalty matrices, and the ordered predictor names through the basis
signatures. Covariance inverses and discrepancy precisions are spectrally
stabilized when required. Nonconvergent path points are retained for
diagnostic purposes but are assigned infinite validation loss and cannot be
selected.
}
\value{
An object of class \code{"cvs_fofr_fit"}. It is a list containing:
  \item{target_fit}{The target-site object supplied in \code{target_fit}.}
  \item{source_summaries}{The validated source summaries used in transfer.}
  \item{gamma_local}{The Target-only coefficient vector of length \eqn{d}.}
  \item{gamma_cvs}{The feasible CVS coefficient vector of length \eqn{d}.}
  \item{gamma_pcvs}{
The selected pCVS coefficient vector of length \eqn{d}, or \code{NULL} when
no pCVS penalty has been selected.
  }
  \item{U_star}{
The \eqn{d \times Kd} control-variate weight matrix in the original basis
coordinates.
  }
  \item{U_star_metric}{
The \eqn{d \times Kd} control-variate weight matrix in Gram-metric
coordinates.
  }
  \item{delta_hat}{
The stacked observed target--source coefficient differences in the original
basis coordinates.
  }
  \item{expected_delta}{
The corresponding stacked differences between the exported plug-in means.
  }
  \item{precision_delta}{
The stabilized precision matrix of the stacked differences in the original
basis coordinates.
  }
  \item{pooled_variance}{
The precision-weighted aggregate covariance matrix in the original basis
coordinates.
  }
  \item{functional_metric}{
The \eqn{d \times d} Gram matrix defining integrated coefficient distance.
  }
  \item{metric_root}{A symmetric square root of \code{functional_metric}.}
  \item{metric_inverse_root}{
The corresponding spectrally stabilized inverse square root.
  }
  \item{delta_hat_metric}{
The stacked observed target--source differences in Gram-metric coordinates.
  }
  \item{expected_delta_metric}{
The stacked plug-in mean differences in Gram-metric coordinates.
  }
  \item{precision_delta_metric}{
The stabilized discrepancy precision used by the pCVS optimization.
  }
  \item{pooled_variance_metric}{
The precision-weighted aggregate covariance matrix in Gram-metric
coordinates.
  }
  \item{zeta_grid}{The decreasing sequence of pCVS penalties that was fitted.}
  \item{delta_path}{
A matrix whose rows are the estimated discrepancy vectors along the pCVS path,
expressed in the original basis coordinates.
  }
  \item{delta_path_metric}{
The corresponding discrepancy path in Gram-metric coordinates.
  }
  \item{pcvs_coefficient_path}{
A matrix whose rows contain the pCVS coefficient vectors in original basis
coordinates, in the order of \code{zeta_grid}.
  }
  \item{validation_mse}{
Integrated validation curve MSE at every path point; a vector of
\code{NA_real_} values when validation data were not supplied.
  }
  \item{selected_zeta}{The selected pCVS penalty, or \code{NULL}.}
  \item{selected_index}{
The index of \code{selected_zeta} in \code{zeta_grid}, or \code{NULL}.
  }
  \item{selection}{A character string describing the selection rule.}
  \item{convergence}{
A logical vector recording convergence of every pCVS path solution.
  }
  \item{iterations}{The number of iterations used at every path point.}
  \item{kkt_residual}{The absolute KKT residual at every path point.}
  \item{relative_kkt_residual}{
The scale-adjusted KKT residual at every path point.
  }
  \item{K}{The number of source studies.}
  \item{group_size}{
The coefficient dimension \eqn{d=(P M_x+Q)M_y} associated with each source
group.
  }
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

x_basis <- fofr_basis("bspline", M = 5)
y_basis <- fofr_basis("bspline", M = 5)

target_fit <- fit_local_fofr(
  target,
  x_basis = x_basis,
  y_basis = y_basis,
  variance_method = "model"
)
source_fit <- fit_local_fofr(
  source,
  x_basis = x_basis,
  y_basis = y_basis,
  variance_method = "model"
)

# Only the aggregate source summary is needed by the target site.
source_summary <- fofr_site_summary(source_fit)

fit <- combine_cvs_fofr(
  target_fit = target_fit,
  source_summaries = list(source_summary),
  target_validation = validation
)

fit$selected_zeta
fit$convergence[fit$selected_index]
coef(fit, method = "CVS")[1:4, 1:3]
coef(fit, method = "pCVS")[1:4, 1:3]
}
}
