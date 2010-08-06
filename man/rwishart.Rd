\name{rwishart}
\alias{rwishart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random deviates from a Wishart distribution }
\description{
  Draws random deviates from a Wishart pdf.
  

}
\usage{
rwishart(N, df, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ Number of random deviates to draw. }
  \item{df}{ Degrees of freedom for Wishart distribution }
  \item{Sigma}{ Mean of the Wishart from which to draw the deviates }
}
\details{
  Draws N matrices of draws from a Wishart with mean Sigma.  This is
  used to draw error covariances for the VAR and BVAR models which are
  distributed inverse Wisharts deviates.
}
\value{

  Returns an N dimensional array of dim(Sigma) square matrices for the
  Wishart random deviates.  If N=1, it returns a single matrix.

}
%\references{ }
\author{ Patrick T. Brandt }
%\note{ }

\seealso{ See also as \code{\link{rmultnorm}} for multivariate normal
  deviates, \code{\link{rgamma}} for the univariate analog to drawing
  Wishart deviates, and \code{\link{ldwishart}} for computing the log
  density for a Wishart variate.}

\examples{
x <- matrix(rnorm(100), 50, 2)
XX <- crossprod(x)
tmp <- rwishart(1, 50, XX)
}
\keyword{ distribution}

