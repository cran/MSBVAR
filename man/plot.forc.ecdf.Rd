\name{plot.forc.ecdf}
\alias{plot.forc.ecdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots VAR forecasts and their empirical error bands}
\description{
  Plots mean VAR forecasts and pointwise error bands 
}
\usage{
\method{plot}{forc.ecdf}(x, probs = c(0.05, 0.95),
               xlab = "", ylab = "", ylim = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ N x nstep matrix of forecasts }
  \item{probs}{ width of error band probabilities, default is 90\% quantiles or
    \code{c(0.05,0.95)}}
  \item{xlab}{ x-axis labels }						  
  \item{ylab}{ y-axis labels }						  
  \item{ylim}{ Bounds for y-axis in standard format \code{c(lower,upper)}}
  \item{\dots}{ other plot  parameters}
  }
\details{
  Plots the mean forecast and the pointwise empirical confidence region
  for a posterior sample of VAR forecasts.
}
\value{
  None.
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ \code{\link{plot.forecast}}}
\examples{
\dontrun{
data(IsraelPalestineConflict)

# Fit a BVAR model
fit.BVAR <- szbvar(IsraelPalestineConflict, p=6, z=NULL, lambda0=0.6,
                   lambda1=0.1, lambda3=2, lambda4=0.5, lambda5=0,
                   mu5=0, mu6=0, nu=3, qm=4, prior=0,
                   posterior.fit=FALSE)

# Generate unconditional forecasts for both models
forecast.BVAR <- uc.forecast(fit.BVAR, nsteps=12,
                                 burnin=100, gibbs=1000)

# Plot the forecasts
par(mfrow=c(2,1))

plot(forecast.BVAR$forecast[,,1], probs=c(0.16,0.84),
               main="I2P Forecast")
abline(h=0)

plot(forecast.BVAR$forecast[,,2], probs=c(0.16,0.84),
               main="P2I Forecast")
abline(h=0)
}
}
\keyword{ hplot}

