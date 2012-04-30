\name{plot.irf}
\alias{plot.irf}
%\alias{plot.irf.VAR}
%\alias{plot.irf.BVAR}
%\alias{plot.irf.BSVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots impulse responses }
\description{
  Plots the m x m matrix of impulse responses produced by \code{irf}.
  
}
\usage{
\method{plot}{irf}(x, varnames = attr(x, "eqnames"), ...)

}

\arguments{
  \item{x}{ Impulse response object produced by \code{\link{irf}}}
  \item{varnames}{ Names of equations and shocks in the format
    \code{c("name1","name2",...)}. Default is to use the names of the
  input variables from the estimation method.}
  \item{...}{other plot arguments}
}
%Improve details to reflect VAR/BVAR/BSVAR availability
\details{
  Generates a plot in ther current plotting device of the impulse
  responses in \code{irf}.  See below for functions that allow one to
  add error bands and confidence regions to the impulse responses.
  Impulses or shocks are in the columns and the rows are the responses.
}
\value{
  None.  Draws a graph in the current device.

}
\references{ Hamilton, James. 1994.  Time Series Analysis, Chapter 11.
  
  Sims, C.A. 1980.  "Macroeconomics and Reality" Econometrica. 
  }
\author{ Patrick T. Brandt}
\note{ This function should NOT be used for Monte Carlo samples of
  IRFs.  Use \code{\link{plot.mc.irf}} for this purpose.}

\seealso{ \code{\link{irf}} to produce impulse responses from a VAR
  object, \code{\link{mc.irf}}, and \code{\link{plot.mc.irf}}
  for methods that allow frequentist and Bayesian error bands in the
  impulse responses}
\examples{
data(IsraelPalestineConflict)
rf.var <- reduced.form.var(IsraelPalestineConflict, p=6)
plot(irf(rf.var, nsteps = 12))
}
\keyword{ hplot}
\keyword{ models}

