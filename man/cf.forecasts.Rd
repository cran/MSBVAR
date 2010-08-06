\name{cf.forecasts}
\alias{cf.forecasts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare VAR forecasts to each other or real data } 
\description{
  Computes the root mean sqaured error and mean absolute error for a
  series of forecasts or for forecasts and real data.  
}
\usage{
cf.forecasts(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{ Matrix of VAR forecasts  produced by \code{forecast.VAR}.}
  \item{m2}{ Matrix of VAR forecasts or a matrix of real data to compare
    to forecasts.}
}
\details{
  Simple RMSE and MAE computation for the forecasts.  The reported values
  are summed over the series and time points.
}
\value{
  An object with two elements:
  \item{rmse }{Forecast RMSE}
  \item{mae }{Forecast MAE}
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{}

\seealso{ \code{\link{forecast}} for forecast computations}
\examples{
data(IsraelPalestineConflict)
Y.sample1 <- window(IsraelPalestineConflict, end=c(2002, 52))
Y.sample2 <- window(IsraelPalestineConflict, start=c(2003,1))

# Fit a BVAR model
fit.bvar <- szbvar(Y.sample1, p=6, lambda0=0.6, lambda1=0.1, lambda3=2,
                   lambda4=0.25, lambda5=0, mu5=0, mu6=0, prior=0)

# Forecast -- this gives back the sample PLUS the forecasts!

forecasts <- forecast(fit.bvar, nsteps=nrow(Y.sample2))

# Compare forecasts to real data
cf.forecasts(forecasts[(nrow(Y.sample1)+1):nrow(forecasts),], Y.sample2)
}
\keyword{ ts }% at least one, from doc/KEYWORDS

