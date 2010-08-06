\name{normalize.svar}
\alias{normalize.svar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Likelihood normalization of SVAR models}
\description{
  Computes various sign normalizations of Bayesian structural VAR (B-SVAR) models.
}
\usage{
normalize.svar(A0unnormalized, A0mode,
               method = c("DistanceMLA", "DistanceMLAhat",
                          "Euclidean", "PositiveDiagA",
                          "PositiveDiagAinv", "Unnormalized"),
               switch.count = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A0unnormalized}{ \eqn{m \times m}{m x m} unnormalized matrix 
    value of \eqn{A_0}{A(0)} in an B-SVAR  }
  \item{A0mode}{ \eqn{m \times m}{m x m} matrix of the \eqn{A_0}{A(0)}
    to normalize around }
  \item{method}{ string that selects the normalization method }
  \item{switch.count}{ counter that counts the number of sign switches.
  Can be non-zero if you want to track the sign switches iteratively.}
}
\details{
  The likelihood of VAR models are invariant to sign changes of the
  structural equation coefficients across equations.  Thus a VAR with
  \eqn{m} equations has a likelihood with \eqn{2^m} identical peaks,
  each a different set of signs (but with the same posterior peak).
  Normalization is used to choose among these peaks.  The most common
  choice is to select the peak where the diagonal elements of
  \eqn{A_0}{A(0)} are all positive, but will not be possible in all
  cases since no such normalization may exist.  Thus, one should select
  a single peak and map all of the draws back to that peak.

  The available normalization methods are
  1) "DistanceMLA" : normalize around the ML peak of A0mode,
  2) "DistanceMLAhat" : normalize around the ML peak of inv(A0mode) 
  3) "Euclidean" : normalize by minimizing the distance between the two matrices.
  4) "PositiveDiagA" : normalize by making the diagonal positive
  5) "PositiveDiagAinv" : normalize by making the diagonal of inv(A0) positive.
  6) "Unnormalized" : no normalization is performed and the function
  returns A0 unnormalized.
}
\value{
  A list with two elements
  \item{A0normalized }{ \eqn{m \times m}{m x m} matrix, the normalized
    value of \eqn{A_0}{A(0)} according to the selected normalization rule.}
  \item{switch.count }{Number of signs changed in the normalization}

}
\references{
  Waggoner, Daniel F. and Tao A. Zha.  2003a. "A Gibbs sampler for
  structural vector autoregressions" \emph{Journal of Economic Dynamics
  \& Control}. 28:349--366.

Waggoner, Daniel F. and Tao A. Zha. 2003b. "Likelihood preserving
normalization in multiple equation models".  \emph{Journal of
Econometrics}. 114: 329--347.
}
\author{ Patrick T. Brandt}
\note{ This function is called in \code{gibbs.A0.BSVAR}, the
  Gibbs sampling of
  \code{szbsvar} models.  In those
  functions, the \eqn{A_0}{A(0)} produced by \code{szbsvar} is
  unnormalized.  The Gibbs sampled draws are
  then normalized using the "DistanceMLA" method, which is consistent
  with the positive system shocks typically seen in the literature, if
  such a normalization exists.  Note that Waggoner and Zha prefer the
  "DistanceMLA" method.
}

\seealso{ \code{\link{szbsvar}}, \code{\link{gibbs.A0}} }
%\examples{}
\keyword{ ts }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
