\name{decay.spec}
\alias{decay.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lag decay specification check} 
\description{ Provides a quick way to visualize the lag decay
  specification in a BVAR model for given parameters by computing the
  variance of the prior VAR coefficients across various lags.  
}
\usage{
decay.spec(qm, p, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qm}{ Periodicity parameter: either 4 or 12 for quarterly or
    monthly data.}
  \item{p}{ Number of lags}
  \item{lambda}{ Lag decay parameter [>0], which is \code{lambda3} in
    the Sims-Zha BVAR specification in \code{szbvar}}
}
\details{
  Computes the relative decay in the prior variance of the VAR prior
  across the lags from 1 to
  p.  Useful for visualizing the rate of decay or how tight the prior
  becomes at higher order lags.
}
\value{
  A time series of length \code{p} of the prior variances for each lag.  

}
\references{ Sims, C.A. and Tao Zha. 1998. "Bayesian Methods for Dynamic
  Multivariate Models." \emph{International Economic Review}. 39(4):949-968.}
\author{ Patrick T. Brandt}
%\note{}

\seealso{ \code{\link{szbvar}}}
\examples{
# Harmonic lag decay example
harmonic <- decay.spec(4, 6, 1)

# Quadratic lag decay example
quadratic <- decay.spec(4, 6, 2)

plot(cbind(harmonic,quadratic))

}
\keyword{ ts}
