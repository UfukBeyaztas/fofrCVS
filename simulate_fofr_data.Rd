\name{simulate_fofr_data}
\alias{simulate_fofr_data}
\title{Simulate Mixed-Predictor Function-on-Function Data}
\description{
Generates a functional response, two noisy functional predictors, and two
scalar covariates from the mixed-predictor function-on-function regression
model used in the package examples. Common and shifted regression-effect
regimes are available for illustrating transferable and nontransferable
source studies.
}
\usage{
simulate_fofr_data(
  n,
  s_grid = seq(0, 1, length.out = 101),
  t_grid = seq(0, 1, length.out = 101),
  covariance_scale = 4,
  effect = c("common", "shifted"),
  measurement_sd = 0.05,
  error_sd = 0.10
)
}
\arguments{
  \item{n}{
A single integer of at least 5 giving the number of subjects.
  }
  \item{s_grid}{
A numeric vector containing at least two finite, strictly increasing
predictor-domain locations.
  }
  \item{t_grid}{
A numeric vector containing at least two finite, strictly increasing
response-domain locations.
  }
  \item{covariance_scale}{
A positive finite number multiplying the covariance kernel of both latent
functional predictors. Changing this value alters the predictor distribution
without changing the regression effects.
  }
  \item{effect}{
The regression-effect regime. \code{"common"} uses the reference coefficient
functions; \code{"shifted"} changes both function-on-function surfaces and
both scalar-effect functions. Partial matching is supported.
  }
  \item{measurement_sd}{
A nonnegative finite number giving the standard deviation of independent
measurement error added to each latent functional-predictor value.
  }
  \item{error_sd}{
A nonnegative finite number giving the standard deviation of independent
response error at each response-grid location.
  }
}
\details{
For each subject \eqn{i}, two independent latent predictor curves,
\eqn{X_{iA}^{\circ}} and \eqn{X_{iE}^{\circ}}, are generated as mean-zero
Gaussian processes. On \code{s_grid}, both have covariance kernel
\deqn{
C(s,u)=a\exp\{-10|s-u|\},
}
where \eqn{a} is \code{covariance_scale}. The implementation adds
\eqn{10^{-9}} to the diagonal of the discretized covariance matrix before its
Cholesky factorization. The two scalar covariates are generated independently
as
\deqn{
w_{i1}\sim N(0,1), \qquad w_{i2}\sim \mathrm{Bernoulli}(0.5).
}
The observed functional predictors are
\deqn{
Z_{ip}(s_j)=X_{ip}^{\circ}(s_j)+e_{ipj},
\qquad e_{ipj}\sim N(0,\sigma_m^2),
}
for \eqn{p\in\{A,E\}}, where \eqn{\sigma_m} is \code{measurement_sd}.
Measurement errors are independent across subjects, predictors, and grid
locations and are independent of the latent curves and scalar covariates.

The functional response is generated from
\deqn{
Y_i(t)=
\int_{\mathcal{S}}X_{iA}^{\circ}(s)\beta_A(s,t)\,ds
+\int_{\mathcal{S}}X_{iE}^{\circ}(s)\beta_E(s,t)\,ds
+w_{i1}\alpha_1(t)+w_{i2}\alpha_2(t)+\varepsilon_i(t).
}
The integrals are evaluated numerically on \code{s_grid} using trapezoidal
weights. The response errors are independent
\eqn{N(0,\sigma_\varepsilon^2)} variables across subjects and response-grid
locations, with \eqn{\sigma_\varepsilon} equal to \code{error_sd}. The
returned \code{conditional_mean} is the noise-free mean conditional on the
latent predictor curves and scalar covariates, rather than on their
measurement-error-contaminated observations.

Under \code{effect = "common"}, the coefficient functions are
\deqn{
\beta_A(s,t)=\sin(2\pi s)\cos(\pi t), \qquad
\beta_E(s,t)=2(s-0.5)(t-0.5),
}
and
\deqn{
\alpha_1(t)=1+0.5\sin(2\pi t), \qquad
\alpha_2(t)=-0.8t.
}
Under \code{effect = "shifted"}, they become
\deqn{
\beta_A^{\mathrm{shift}}(s,t)=-0.35\beta_A(s,t),
}
\deqn{
\beta_E^{\mathrm{shift}}(s,t)=\beta_E(s,t)
+1.25\cos(\pi s)\sin(2\pi t),
}
and
\deqn{
\alpha_1^{\mathrm{shift}}(t)=-0.50\alpha_1(t), \qquad
\alpha_2^{\mathrm{shift}}(t)=-0.8t+0.9\sin(\pi t).
}
Thus, independently generated studies with \code{effect = "common"} share
the same regression operator and are transferable under the simulation model.
Changing only \code{covariance_scale} creates predictor-covariance shift while
preserving that operator. By contrast, \code{effect = "shifted"} creates a
source with a genuinely different regression operator and is intended for
illustrating possible negative transfer and the protective role of pCVS.
}
\value{
A list containing:
  \item{Y}{
An \eqn{n\times J_t} matrix of observed response curves, where \eqn{J_t} is
the length of \code{t_grid}.
  }
  \item{y_grid}{The response grid supplied in \code{t_grid}.}
  \item{X_list}{
A named list with elements \code{activity} and \code{exposure}. Each is an
\eqn{n\times J_s} matrix of measurement-error-contaminated functional
predictors, where \eqn{J_s} is the length of \code{s_grid}.
  }
  \item{x_grids}{
A named list assigning \code{s_grid} to both functional predictors.
  }
  \item{w}{
An \eqn{n\times 2} matrix with columns \code{continuous} and \code{binary}.
  }
  \item{conditional_mean}{
An \eqn{n\times J_t} matrix containing the noise-free response means
conditional on the latent functional predictors and scalar covariates.
  }
  \item{latent_X_list}{
A named list containing the latent, measurement-error-free \code{activity}
and \code{exposure} predictor matrices.
  }
  \item{truth}{
A nested list. \code{truth$beta} contains the \eqn{J_s\times J_t} matrices
\code{activity} and \code{exposure}; \code{truth$alpha} contains the
length-\eqn{J_t} vectors \code{continuous} and \code{binary}.
  }
  \item{effect}{The selected regression-effect regime.}
  \item{covariance_scale}{The predictor-covariance scale used to generate the
latent curves.}
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_grid <- seq(0, 1, length.out = 101)
t_grid <- seq(0, 1, length.out = 101)

target <- simulate_fofr_data(
  40,
  s_grid = s_grid,
  t_grid = t_grid
)
transferable_source <- simulate_fofr_data(
  80,
  s_grid = s_grid,
  t_grid = t_grid,
  covariance_scale = 7
)
shifted_source <- simulate_fofr_data(
  80,
  s_grid = s_grid,
  t_grid = t_grid,
  effect = "shifted"
)
validation <- simulate_fofr_data(
  40,
  s_grid = s_grid,
  t_grid = t_grid
)

dim(target$Y)
names(target$X_list)
names(target$truth$beta)
mean((target$Y - target$conditional_mean)^2)

fit <- fit_cvs_fofr(
  datasets = list(target, transferable_source, shifted_source),
  x_basis_type = "bspline",
  y_basis_type = "bspline",
  Mx = 5,
  My = 5,
  variance_method = "model",
  target_validation = validation
)

fit
}
}
