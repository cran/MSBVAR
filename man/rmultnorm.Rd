\name{rmultnorm}
\alias{rmultnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate Normal Random Number Generator }
\description{
  Generates multivariate normal random variates for give mean and
  covariance vectors.  Can also handle generation of multivariate normal
  deviates with singular covariance distributions via singular value
  decomposition (SVD).
}
\usage{
rmultnorm(n, mu, vmat, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Number of variates to draw.}
  \item{mu}{ M column matrix of multivariate means}
  \item{vmat}{ M x M covariance matrix }
  \item{tol}{ Tolerance level used for SVD of the covariance.  Default is 1e-10}
}
\details{
  Generates N draws from a multivariate normal distribution with mean
  matrix mu and covariance matrix vmat.
}
\value{
  Matrix of the random draw that is conformable with the input mu.
}
\references{ }
\author{ Patrick T. Brandt }
\note{ Based on code by Jeff Gill. This function is called in the hard
  condition forecasting in \code{\link{hc.forecast.var}} for simulating
  the structural innovations.}

\seealso{ \code{\link{rnorm}} }
\examples{
rmultnorm(1, matrix(c(1,2),2,1), vmat=matrix(c(1,1,0,1),2,2))
}
\keyword{ distribution}

