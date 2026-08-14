\name{coefficient_surface}
\alias{coefficient_surface}
\title{Evaluate a Fitted Function-on-Function Coefficient Surface}
\description{
Reconstructs a bivariate coefficient surface for one functional predictor
from the Target-only, control-variates (CVS), or penalized control-variates
(pCVS) estimator in a fitted function-on-function transfer model.
}
\usage{
coefficient_surface(
  object,
  predictor = 1,
  method = c("local", "CVS", "pCVS"),
  s_grid = NULL,
  t_grid = NULL
)
}
\arguments{
  \item{object}{
An object of class \code{"cvs_fofr_fit"}, as returned by
\code{\link{fit_cvs_fofr}} or \code{\link{combine_cvs_fofr}}.
  }
  \item{predictor}{
The functional predictor whose coefficient surface is to be evaluated. Supply
either a single integer between 1 and the number of functional predictors or
the exact predictor name recorded in \code{object}.
  }
  \item{method}{
A character string specifying the estimator. \code{"local"} uses the
Target-only estimator, \code{"CVS"} uses the feasible control-variates
estimator, and \code{"pCVS"} uses the selected penalized control-variates
estimator. Partial matching is supported.
  }
  \item{s_grid}{
A numeric vector of predictor-domain locations at which the surface is to be
evaluated. If \code{NULL}, 101 equally spaced points over the boundary of the
fitted predictor basis are used.
  }
  \item{t_grid}{
A numeric vector of response-domain locations at which the surface is to be
evaluated. If \code{NULL}, 101 equally spaced points over the boundary of the
fitted response basis are used.
  }
}
\details{
For the \eqn{p}th functional predictor, the function-on-function linear model
contains the term
\deqn{
\int_{\mathcal{S}} X_{ip}(s)\,\beta_p(s,t)\,ds,
}
so \eqn{\beta_p(s,t)} describes how the predictor value at location \eqn{s}
contributes to the response at location \eqn{t}. The fitted model represents
this surface in a tensor-product sieve. If
\eqn{\boldsymbol{\phi}(s)} is the predictor basis,
\eqn{\boldsymbol{\vartheta}(t)} is the response basis, and
\eqn{\widehat{\boldsymbol{\Gamma}}_p} is the coefficient block associated
with predictor \eqn{p}, the reconstructed surface is
\deqn{
\widehat{\beta}_p(s,t)
=\boldsymbol{\phi}(s)^\mathsf{T}
 \widehat{\boldsymbol{\Gamma}}_p
 \boldsymbol{\vartheta}(t).
}

The function first selects the requested Target-only, CVS, or pCVS coefficient
vector, reshapes it into the fitted coefficient matrix, and extracts the
\eqn{M_x \times M_y} block belonging to \code{predictor}. It then evaluates
the two bases on \code{s_grid} and \code{t_grid} and performs the
tensor-product reconstruction above.

The supplied grids are evaluation grids only. They may differ from the grids
on which the functional observations were recorded, and changing them does not
refit the model. For the usual interpretation, their values should lie within
the boundaries of the fitted predictor and response bases. The returned
surface has one row for each value of \code{s_grid} and one column for each
value of \code{t_grid}; explicitly,
\code{beta[i, j]} is the estimate at
\code{(s_grid[i], t_grid[j])}.

Although the numerical basis coefficients depend on the chosen bases, the
reconstructed coefficient surface is the functional object represented by
those coefficients. In particular, the Gram-metric formulation used by pCVS
ensures that its discrepancy penalty measures integrated surface distance
rather than an arbitrary Euclidean distance between basis coefficients.
Extraction with \code{method = "pCVS"} requires a selected pCVS solution;
otherwise an error is returned.
}
\value{
A list with the following components:
  \item{s_grid}{The predictor-domain evaluation grid.}
  \item{t_grid}{The response-domain evaluation grid.}
  \item{beta}{
A numeric matrix of dimension \code{length(s_grid)} by
\code{length(t_grid)}. Entry \code{beta[i, j]} equals the reconstructed
coefficient surface at \code{(s_grid[i], t_grid[j])}.
  }
  \item{method}{The estimator used for reconstruction.}
  \item{predictor}{The numerical index of the selected functional predictor.}
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

surface <- coefficient_surface(
  fit,
  predictor = "activity",
  method = "pCVS"
)

dim(surface$beta)
surface$beta[1:3, 1:3]

persp(
  surface$s_grid,
  surface$t_grid,
  surface$beta,
  theta = 35,
  phi = 25,
  expand = 0.6,
  shade = 0.4,
  col = "lightblue",
  ticktype = "detailed",
  xlab = "s",
  ylab = "t",
  zlab = expression(hat(beta)(s, t))
)
}
}
