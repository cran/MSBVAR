\name{rmse}
\alias{rmse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Root mean squared error of a Monte Carlo / MCMC sample of forecasts}
\description{
  Computes the root mean squared error (RMSE) of a Monte Carlo sample of
  forecasts. 
}
\usage{
rmse(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{ Forecast sample for model 1}
  \item{m2}{ Forecast sample for model 2}
}
\details{
  User needs to subset the forecasts if necessary.
}
\value{
  Forecast RMSE.
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ \code{\link{mae}}, \code{\link{forecast}}}
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
rmse(forecasts[(nrow(Y.sample1)+1):nrow(forecasts),], Y.sample2)

}
\keyword{ ts}

