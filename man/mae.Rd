\name{mae}
\alias{mae}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mean absolute error of VAR forecasts }
\description{
  Computes the mean absolute error of VAR forecasts
}
\usage{
mae(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{ \eqn{nsteps \times m}{nstep x m} matrix of VAR forecasts }
  \item{m2}{ \eqn{nsteps \times m}{nstep x m} matrix of VAR forecasts or true values}
}
\details{
  Computes the mean absolute error (MAE) across a series of VAR forecasts.
}
\value{
  MAE value
}
\author{ Patrick T. Brandt}
%\note{ }
\seealso{ \code{\link{cf.forecasts}}, \code{\link{rmse}}}
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
mae(forecasts[(nrow(Y.sample1)+1):nrow(forecasts),], Y.sample2)

}
\keyword{ ts}
