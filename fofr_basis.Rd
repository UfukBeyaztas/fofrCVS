\name{fofr_basis}
\alias{fofr_basis}
\title{Construct a Basis for Function-on-Function Regression}
\description{
Constructs a normalized Fourier or B-spline basis together with its numerical
Gram matrix and integrated second-derivative penalty matrix. The resulting
object supplies the common functional geometry used by local FoFR fitting,
CVS aggregation, and Gram-metric pCVS regularization.
}
\usage{
fofr_basis(
  type = c("fourier", "bspline"),
  M = 15,
  degree = 3,
  boundary = c(0, 1),
  integration_grid = seq(boundary[1], boundary[2], length.out = 1001)
)
}
\arguments{
  \item{type}{
The basis family: \code{"fourier"} or \code{"bspline"}. Partial matching is
supported.
  }
  \item{M}{
The number of basis functions. It must be a single integer of at least 3. For
a B-spline basis, it must also satisfy \code{M >= degree + 1}.
  }
  \item{degree}{
The polynomial degree of a B-spline basis. It must be a single integer of at
least 2. The value is stored but does not affect a Fourier basis.
  }
  \item{boundary}{
A length-two increasing numeric vector defining the functional domain.
  }
  \item{integration_grid}{
A strictly increasing numeric vector of at least three finite points lying
within \code{boundary}. It is used for trapezoidal approximation of the Gram
and roughness matrices. For a B-spline basis, it is also the reference grid
used to place the interior knots.
  }
}
\details{
Let \eqn{\boldsymbol{\phi}(u)=
(\phi_1(u),\ldots,\phi_M(u))^\mathsf{T}} denote the constructed basis on
\eqn{[a,b]}, where \eqn{a} and \eqn{b} are the entries of
\code{boundary}. The function returns numerical approximations to
\deqn{
\boldsymbol{\Psi}
=\int_a^b\boldsymbol{\phi}(u)\boldsymbol{\phi}(u)^\mathsf{T}\,du,
\qquad
\boldsymbol{W}
=\int_a^b\boldsymbol{\phi}^{\prime\prime}(u)
\boldsymbol{\phi}^{\prime\prime}(u)^\mathsf{T}\,du.
}
For a coefficient vector \eqn{\boldsymbol{c}} representing
\eqn{f(u)=\boldsymbol{\phi}(u)^\mathsf{T}\boldsymbol{c}}, these matrices
satisfy
\deqn{
\int_a^b f(u)^2\,du
=\boldsymbol{c}^\mathsf{T}\boldsymbol{\Psi}\boldsymbol{c},
\qquad
\int_a^b\{f^{\prime\prime}(u)\}^2\,du
=\boldsymbol{c}^\mathsf{T}\boldsymbol{W}\boldsymbol{c}.
}
Thus, \code{Psi} represents the \eqn{L^2} geometry of the basis and
\code{W} represents second-derivative roughness. For a tensor-product
coefficient surface
\eqn{\beta(s,t)=\boldsymbol{\phi}(s)^\mathsf{T}
\boldsymbol{\Gamma}\boldsymbol{\vartheta}(t)}, the corresponding integrated
squared norm is
\deqn{
\int\!\!\int \beta(s,t)^2\,ds\,dt
=\mathrm{tr}(\boldsymbol{\Gamma}^\mathsf{T}
\boldsymbol{\Psi}_x\boldsymbol{\Gamma}\boldsymbol{\Psi}_y).
}
This identity is the basis of the Gram metric used to compare entire
coefficient objects in pCVS.

For \code{type = "fourier"}, put \eqn{D=b-a} and
\eqn{v=(u-a)/D}. The basis is ordered as
\deqn{
\phi_1(u)=D^{-1/2},\qquad
\phi_{2r}(u)=(2/D)^{1/2}\cos(2\pi r v),\qquad
\phi_{2r+1}(u)=(2/D)^{1/2}\sin(2\pi r v),
}
for the frequencies represented among the first \eqn{M} functions. This
normalization makes the continuous Gram matrix the identity. The constant
function has zero second-derivative penalty, while the trigonometric penalties
increase with the fourth power of frequency.

For \code{type = "bspline"}, the function constructs an intercept-containing
B-spline basis with \code{M} functions and the requested polynomial degree.
Interior knots are those selected by \code{splines::bs} from
\code{integration_grid}; the full knot sequence repeats each boundary
\code{degree + 1} times. The same knot sequence is subsequently used to
evaluate the basis and its second derivatives on study-specific observation
grids.

Both \code{Psi} and \code{W} are computed by trapezoidal quadrature on
\code{integration_grid} and are explicitly symmetrized. Increasing the grid
density generally improves their approximation to the continuous integrals.
For Fourier bases, the basis functions themselves do not depend on this grid.
For B-spline bases, changing a nonuniform integration grid may also change the
automatically placed interior knots and therefore changes the basis itself.
}
\value{
An object of class \code{"fofr_basis"}. It is a list containing:
  \item{type}{The selected basis family.}
  \item{M}{The number of basis functions.}
  \item{degree}{The supplied polynomial degree.}
  \item{boundary}{The two endpoints of the functional domain.}
  \item{integration_grid}{The grid used for knot placement and quadrature.}
  \item{knots}{
For a B-spline basis, the interior knots returned by \code{splines::bs}. This
component is absent for a Fourier basis.
  }
  \item{full_knots}{
For a B-spline basis, the complete knot sequence used for basis evaluation,
including repeated boundary knots. This component is absent for a Fourier
basis.
  }
  \item{Psi}{
The symmetric \eqn{M\times M} numerical Gram matrix.
  }
  \item{W}{
The symmetric \eqn{M\times M} integrated second-derivative penalty matrix.
  }
}
\author{
Ufuk Beyaztas, Muge Mutis, Sude Gurer and Soutir Bandyopadhyay.
}
\examples{
fourier_basis <- fofr_basis("fourier", M = 7)
bspline_basis <- fofr_basis("bspline", M = 7, degree = 3)

dim(fourier_basis$Psi)
round(fourier_basis$Psi, 4)
round(diag(fourier_basis$W), 3)

bspline_basis$knots
range(eigen(
  bspline_basis$Psi,
  symmetric = TRUE,
  only.values = TRUE
)$values)

stopifnot(
  isSymmetric(fourier_basis$Psi),
  isSymmetric(fourier_basis$W),
  isSymmetric(bspline_basis$Psi),
  isSymmetric(bspline_basis$W)
)
}
