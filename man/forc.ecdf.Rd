\name{forc.ecdf}
\alias{forc.ecdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Empirical CDF computations for posterior forecast samples }
\description{
  Computes (pointwise over time) empirical density (error bands) and
  mean forecasts for a Monte Carlo or Bayesian posterior sample of forecasts.
}
\usage{
forc.ecdf(forecasts, probs = c(0.05, 0.95), start = c(0, 1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forecasts}{ Posterior sample of VAR forecasts produced by
    \code{hc.forecast.VAR()} or \code{uc.forecast.VAR()}}
  \item{probs}{ Error band width in percentiles, default is 90\% error band. }
  \item{start}{ Start value for the time series -- as in the \code{ts()}
  for the forecast horizon}
  \item{\dots}{ Other \code{ecdf()} parameters}
}
\details{
  For each endogenous variable in the VAR and each point in the forecast
  horizon this function estimates the percentile based confidence
  interval.  It then returns a time series matrix beginning at
  \code{start} of the mean forecast and the limits of the confidence
  region for each variable in the forecast sample.
}
\value{
  A multiple time series object is returned where the first column is
  the mean estimate followed by the upper and lower bounds of the
  confidence region.
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{}

%\seealso{ }
%\examples{}
\keyword{ ts }% at least one, from doc/KEYWORDS

