\name{plot.forecast}
\alias{plot.forecast}
\alias{plot.forecast.VAR}
\alias{plot.forecast.BVAR}
\alias{plot.forecast.BSVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots competing sets of VAR forecasts or a single set of VAR forecasts}
\description{
  Produces a high level plot of two sets of VAR forecasts computed from
  a Markov Chain Monte Carlo (MCMC) posterior sample.  Forecasts are the
  output of \code{hc.forecast}, \code{forecast}, and \code{uc.forecast}.
}
\usage{
plot.forecast.VAR(x, y = NULL, varnames = NULL, start = c(0, 1), freq = 1, probs =
     c(0.05, 0.95), compare.level = NULL, ylab = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ First sample of forecasts }
  \item{y}{ Second sample of forecasts, default = \code{NULL}}
  \item{varnames}{ m x 1 list of variable names, e.g.,
    \code{c("name1","name2",\dots)}}
  \item{start}{ Start date for the forecast period time series, default
    = \code{c(0,1)}}
  \item{freq}{ Time series frequency, in format consistent with
    \code{ts}.  Default is 1}
  \item{probs}{ Probability limits for error bands.  Default is 90\% or
    \code{c(0.05,0.95)}}
  \item{compare.level}{ ~~Describe \code{compare.level} here~~ }
  \item{ylab}{ ~~Describe \code{ylab} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  Plots the mean forecast and the pointwise empirical confidence region
  for a posterior sample of VAR forecasts.  Overlays a second set of
  forecasts and error bands if requested in \code{fcasts2}
}
\value{
  None. Plots forecasts on the current display device.
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{}

\seealso{ \code{\link{plot.forc.ecdf}}, and
  \code{\link{uc.forecast}} for an example.}
%\examples{}
\keyword{ hplot}
\keyword{ models}
