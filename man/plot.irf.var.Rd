\name{plot.irf.var}
\alias{plot.irf.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots impulse responses }
\description{
  Plots the m x m matrix of impulse responses produced by \code{irf.var}.
}
\usage{
plot.irf.var(x, varnames = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Impulse response object produced by \code{irf.var}}
  \item{varnames}{ Variance names in the format
    \code{c("name1","name2",...)} }
  \item{...}{other plot arguments}
}
\details{
  Generates a plot in ther current plotting device of the impulse
  responses in \code{irf}.  See below for functions that allow one to
  add error bands and confidence regions to the impulse responses.
  Impulses or shocks are in the columns and the rows are the responses.}  
}
\value{
  None.  Draws a graph in the current device.

}
\references{ Hamilton, James. 1994.  Time Series Analysis, Chapter 11.
  
  Sims, C.A. 1980.  "Macroeconomics and Reality" Econometrica. 
  }
\author{ Patrick T. Brandt}
\note{ This function should NOT be used for Monte Carlo samples of
  IRFs.  Use \code{\link{plot.mc.irf.var}} for this purpose.}

\seealso{ \code{\link{irf.var}} to produce impulse responses from a VAR
  object, \code{\link{mc.irf.var}}, and \code{\link{plot.mc.irf.var}}
  for methods that allow frequentist and Bayesian error bands in the
  impulse responses}
\examples{
data(IsraelPalestineConflict)
rf.var <- reduced.form.var(IsraelPalestineConflict, p=6)
plot(irf.var(rf.var, nsteps = 12))
}
\keyword{ hplot}
\keyword{ models}

