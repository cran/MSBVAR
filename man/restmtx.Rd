\name{restmtx}
\alias{restmtx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utility function for generating the restriction matrix for hard
  condition forecasting}
\description{
  Generates the restriction matrix for a set of hard condition
  forecasts.  See \code{\link{hc.forecast}} for details.
}
\usage{
restmtx(nsteps, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsteps}{ Number of periods in the forecast horizon }
  \item{m}{ Number of endogenous variables in the VAR.}
}
\details{
  Builds the appropriately dimensioned and filled restriction  matrix of
  zeros and ones for hard condition forecasting.
}
\value{
  A matrix of dimensions (nsteps x m*nsteps) that can be used to
  represent the restrictions in hard condition forecasting using
  \code{hc.forecast}
}

\references{ Waggoner, Daniel F. and Tao Zha. 1999. "Conditional
  Forecasts in Dynamic Multivariate Models" \emph{Review of Economics and
  Statistics}, 81(4):639-651.}
\author{ Patrick T. Brandt}

\seealso{ \code{\link{hc.forecast}}}
%\examples{}
\keyword{ models}
\keyword{ utilities}
