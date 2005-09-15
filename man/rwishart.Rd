\name{rwishart}
\alias{rwishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random deviates from a Wishart distribution }
\description{
  Draws deviates from a Wishart pdf.
}
\usage{
rwishart(N, df, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ Number of deviates to draw. }
  \item{df}{ Degrees of reedom for Wishart distribution }
  \item{Sigma}{ Mean of the Wishart from which to draw the deviates }
}
\details{
  Draws N matrices of drawss from a Wishart with mean Sigma.  
}
\value{

  Returns an N dimensional array of dim(Sigma) matrices for the Wishart random
  deviates.  If N=1, it returns a single matrix.

}
\references{ }
\author{ Patrick T. Brandt }
\note{ }

\seealso{ See also as \code{\link{rmultnorm}} for multivariate normal
  deviates and \code{\link{rgamma}} for the univariate analog to drawing
  Wishart deviates.}
\examples{
x <- matrix(rnorm(100), 50, 2)
XX <- crossprod(x)
tmp <- rwishart(1, 50, XX)

}
\keyword{ distribution}

