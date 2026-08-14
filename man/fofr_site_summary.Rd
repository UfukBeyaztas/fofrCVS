\name{fofr_site_summary}
\alias{fofr_site_summary}
\title{Create a Transfer Summary from a Local FoFR Fit}
\description{
Extracts the aggregate coefficient, moment, and basis information required to
use a locally fitted function-on-function regression model as a source in the
CVS and pCVS transfer procedures. Subject-level responses, predictors, and
covariates are not included in the returned object.
}
\usage{
fofr_site_summary(local_fit)
}
\arguments{
  \item{local_fit}{
An object of class \code{"local_fofr_fit"} produced by
\code{\link{fit_local_fofr}}.
  }
}
\details{
The function defines the communication boundary between a source site and the
coordinating target site. Let \eqn{P} be the number of functional predictors,
\eqn{Q} the number of scalar covariates, and \eqn{M_x} and \eqn{M_y} the
predictor- and response-basis dimensions. The local coefficient matrix has
dimension \eqn{(P M_x+Q)\times M_y}, and the length of its columnwise
vectorization is
\deqn{
d=(P M_x+Q)M_y.
}
For study \eqn{k}, the exported numerical summary consists of
\deqn{
\mathcal{S}_k=
\left(\widehat{\boldsymbol{g}}_k,
\widehat{\boldsymbol{\mu}}_k,
\widehat{\boldsymbol{\Sigma}}_k\right),
}
where
\eqn{\widehat{\boldsymbol{g}}_k=
\mathrm{vec}(\widehat{\boldsymbol{\Gamma}}_k)} is the local coefficient
vector, \eqn{\widehat{\boldsymbol{\mu}}_k} is its feasible plug-in
conditional-mean vector, and \eqn{\widehat{\boldsymbol{\Sigma}}_k} is the
selected plug-in covariance estimate. The mean and covariance are the objects
stored as \code{expected_gamma} and \code{variance_gamma} by
\code{\link{fit_local_fofr}}. In particular, \code{expected_gamma} should not
in general be interpreted as an exact conditional expectation; its exactness
depends on the working model and on how the smoothing parameters are treated.

In addition to these numerical summaries, the returned object records the
basis dimensions and a basis signature. The signature encodes the predictor
and response basis definitions, their Gram and roughness matrices, and the
ordered functional-predictor and scalar-covariate names. These checks are
essential: target--source coefficient subtraction is meaningful only when all
sites use the same coefficient coordinates and the same ordering of model
terms. Observation grids may differ across sites, provided that they lie
within the common basis boundaries and the fitted models share an identical
basis signature.

The numerical communication burden is \eqn{2d+d^2} values---two vectors of
length \eqn{d} and one \eqn{d\times d} covariance matrix---apart from the
small dimension and signature metadata. It is therefore independent of the
number of subjects at the source site. The object does not contain the raw
response matrix, functional-predictor matrices, scalar covariates,
subject-specific residuals, observation grids, smoothing operators, or tuning
paths.

This function performs no cross-study aggregation and does not itself compute
a transfer estimator. A list of source summaries is supplied to
\code{\link{combine_cvs_fofr}}, together with the target-site local fit, to
construct the CVS estimator and the pCVS path.
}
\value{
An object of class \code{"cvs_fofr_summary"}. It is a list containing:
  \item{gamma_hat}{
The local coefficient vector
\eqn{\widehat{\boldsymbol{g}}=
\mathrm{vec}(\widehat{\boldsymbol{\Gamma}})}, of length \eqn{d}.
  }
  \item{expected_gamma}{
The feasible plug-in conditional-mean vector of length \eqn{d} used to center
the target--source coefficient differences.
  }
  \item{variance_gamma}{
The spectrally regularized \eqn{d\times d} plug-in covariance estimate for
\code{gamma_hat}. Its construction is determined by the
\code{variance_method} used in \code{\link{fit_local_fofr}}.
  }
  \item{Mx}{The predictor-basis dimension \eqn{M_x}.}
  \item{My}{The response-basis dimension \eqn{M_y}.}
  \item{P_number}{The number \eqn{P} of functional predictors.}
  \item{Q_number}{The number \eqn{Q} of scalar covariates.}
  \item{basis_signature}{
Metadata identifying the common predictor and response bases and the ordered
model terms. It is used by \code{\link{combine_cvs_fofr}} to reject
incompatible source summaries.
  }
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)
source_data <- simulate_fofr_data(
  60,
  s_grid = s_grid,
  t_grid = t_grid
)

x_basis <- fofr_basis("bspline", M = 5)
y_basis <- fofr_basis("bspline", M = 5)

source_fit <- fit_local_fofr(
  dataset = source_data,
  x_basis = x_basis,
  y_basis = y_basis,
  variance_method = "model"
)

source_summary <- fofr_site_summary(source_fit)

class(source_summary)
names(source_summary)

d <- (source_summary$P_number * source_summary$Mx +
  source_summary$Q_number) * source_summary$My
c(
  coefficient_length = length(source_summary$gamma_hat),
  expected_length = length(source_summary$expected_gamma),
  theoretical_length = d
)
dim(source_summary$variance_gamma)

# No subject-level observations are retained.
intersect(c("Y", "X_list", "w", "x_grids", "y_grid"),
          names(source_summary))
}
}
