\name{predict.cvs_fofr_fit}
\alias{predict.cvs_fofr_fit}
\title{Predict Response Curves from a CVS-FoFR Fit}
\description{
Computes fitted response curves for new functional and scalar predictors from
the Target-only, control-variates (CVS), or penalized control-variates (pCVS)
estimator in a fitted function-on-function transfer model.
}
\usage{
\method{predict}{cvs_fofr_fit}(
  object,
  newdata,
  method = c("local", "CVS", "pCVS"),
  ...
)
}
\arguments{
  \item{object}{
An object of class \code{"cvs_fofr_fit"}, as returned by
\code{\link{fit_cvs_fofr}} or \code{\link{combine_cvs_fofr}}.
  }
  \item{newdata}{
A list containing \code{X_list}, the new functional-predictor matrices.
The list may also contain \code{x_grids}, \code{w}, and \code{y_grid}; see
Details. The scalar-covariate matrix \code{w} is required when the fitted model
contains scalar covariates.
  }
  \item{method}{
A character string specifying the estimator used for prediction.
\code{"local"} uses the Target-only estimator, \code{"CVS"} uses the feasible
control-variates estimator, and \code{"pCVS"} uses the selected penalized
control-variates estimator. Partial matching is supported.
  }
  \item{...}{
Additional arguments for the \code{predict} generic. They are currently
unused.
  }
}
\details{
Suppose the fitted model contains \eqn{P} functional predictors and \eqn{Q}
scalar covariates. For a new subject \eqn{i}, let
\eqn{\widehat{\boldsymbol{h}}_{i0}} be the mixed-predictor design vector formed
using the target-site centering and smoothing quantities. For estimator
\eqn{m\in\{\mathrm{local},\mathrm{CVS},\mathrm{pCVS}\}}, the predicted curve is
\deqn{
\widehat{Y}^{(m)}_i(t)
=\overline{Y}_0(t)
+\widehat{\boldsymbol{h}}_{i0}^{\mathsf{T}}
 \widehat{\boldsymbol{\Gamma}}^{(m)}
 \boldsymbol{\vartheta}(t),
}
where \eqn{\overline{Y}_0(t)} is the target training-sample mean response,
\eqn{\boldsymbol{\vartheta}(t)} is the fitted response basis, and
\eqn{\widehat{\boldsymbol{\Gamma}}^{(m)}} is the coefficient matrix selected
by \code{method}. Thus, transfer modifies the centered regression effects;
all three prediction rules remain anchored by the target-site mean response.

For functional predictor \eqn{p}, the corresponding block of
\eqn{\widehat{\boldsymbol{h}}_{i0}} is
\eqn{\widehat{\boldsymbol{c}}_{ip,0}^{\mathsf{T}}
\boldsymbol{\Psi}_x}, where
\eqn{\widehat{\boldsymbol{c}}_{ip,0}} is obtained after subtracting the target
mean predictor curve and smoothing on the new predictor grid with the
target-fitted basis and smoothing penalty. The scalar block is
\eqn{(\boldsymbol{w}_i-\overline{\boldsymbol{w}}_0)^{\mathsf{T}}}. In
particular, \code{method = "local"} means prediction with the Target-only
coefficient estimator; it does not refit the model to \code{newdata}.

The component \code{X_list} must contain the same number of functional
predictors used to fit \code{object}, with a common number of rows. Each
matrix has one row per new subject and one column per location in the
corresponding grid. A vector is accepted as a single-subject predictor, and a
matrix may be supplied directly when the model has one functional predictor.
If \code{X_list} is fully named, its names must match the fitted predictor
names and its elements are reordered accordingly; otherwise its order is
assumed to match the fitted model. All predictor values must be finite.

The optional \code{x_grids} component supplies one finite, strictly increasing
grid per functional predictor. Its names, when complete, must match the fitted
predictor names. If \code{x_grids} is omitted, the target training grids are
used. New predictor grids may differ from the training grids but must remain
within their respective target training-grid ranges; extrapolation is not
performed. The number of columns of each predictor matrix must equal the
length of its grid.

When \eqn{Q>0}, \code{w} must be a finite matrix with one row per new subject
and \eqn{Q} columns. A vector is reshaped to have the same number of rows as
the functional predictors and must therefore have compatible length. Fully
specified column names must match the fitted scalar covariate names, in which
case the columns are placed in the fitted order. The component is unnecessary
when the model has no scalar covariates.

The optional \code{y_grid} determines the response locations at which the
predictions are evaluated. It must contain at least two finite, strictly
increasing values within the target response training-grid range. If it is
omitted, the target training response grid is used. On a different admissible
grid, the response basis is evaluated anew and the target mean response is
linearly interpolated. A response matrix \code{Y} may be present in
\code{newdata}, for example when a complete test data set is supplied, but it
is not used in forming the predictions.

For \code{method = "CVS"}, the coefficient matrix is obtained from the
feasible control-variates adjustment of the Target-only estimator. For
\code{method = "pCVS"}, it is the selected member of the group-penalized
transfer path. A pCVS prediction is available only after a penalty has been
selected by a target validation sample or supplied through \code{zeta}; an
error is returned otherwise.
}
\value{
A numeric matrix with one row for each new subject and one column for each
location in the selected response grid. If \code{newdata$y_grid} is omitted,
the number of columns equals the length of the target training response grid.
Entry \code{[i, j]} is the predicted response for subject \code{i} at the
\code{j}th response-grid location.
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
\donttest{

s_train <- seq(0, 1, length.out = 101)
t_train <- seq(0, 1, length.out = 101)
s_test <- seq(0, 1, length.out = 101)
t_test <- seq(0, 1, length.out = 101)

target <- simulate_fofr_data(40, s_grid = s_train, t_grid = t_train)
source <- simulate_fofr_data(80, s_grid = s_train, t_grid = t_train)
validation <- simulate_fofr_data(40, s_grid = s_train, t_grid = t_train)
test <- simulate_fofr_data(50, s_grid = s_test, t_grid = t_test)

fit <- fit_cvs_fofr(
  datasets = list(target, source),
  x_basis_type = "bspline",
  y_basis_type = "bspline",
  Mx = 5,
  My = 5,
  variance_method = "model",
  target_validation = validation
)

prediction_data <- list(
  X_list = test$X_list,
  x_grids = test$x_grids,
  w = test$w,
  y_grid = test$y_grid
)

pred_target <- predict(fit, prediction_data, method = "local")
pred_cvs <- predict(fit, prediction_data, method = "CVS")
pred_pcvs <- predict(fit, prediction_data, method = "pCVS")

dim(pred_pcvs)
pred_pcvs[1, 1:5]

# Integrated squared prediction error on a test sample.
response_weights <- c(0.5, rep(1, length(t_test) - 2), 0.5) *
  diff(range(t_test)) / (length(t_test) - 1)
mean(as.numeric((test$Y - pred_pcvs)^2 \%*\% response_weights))
}
}
