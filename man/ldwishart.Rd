\name{ldwishart}
\alias{ldwishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log density for a Wishart variate}

\description{
    Computes log density for a Wishart random variable.
}
\usage{
ldwishart(W, v, S)
}
\arguments{
  \item{W}{ Wishart variate for which the log density is to be computed}
  \item{v}{ degrees of freedom for the Wishart variate}
  \item{S}{ scale factor for the Wishart variate (typically the inverse
    covariance if you are working with a multivariate random normal setup}
}
\details{
  Computes the log density for a Wishart variate with mean \eqn{S} and
    degrees of freedom \eqn{v}.  Special care has been taken to avoid
    underflow in the computation.
}
\value{  A scalar, the value of the log density for the variate \eqn{W}
    with mean \eqn{S} and degrees of freedom \eqn{v}.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{ Patrick T. Brandt }
\note{
  This is modifed from the log density function in MCMCpack.  It better
    handles underflows.
}

\seealso{ \code{\link{rwishart}}}

\examples{
x <- matrix(rnorm(100), 50, 2)
XX <- crossprod(x)
ldwishart(solve(XX), 50, diag(2))
}

\keyword{ distribution }
