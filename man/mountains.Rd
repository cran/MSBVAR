\name{mountains}
\alias{mountains}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mountain plots for summarizing forecast densities }
\description{
  "Mountain plots" summarize the bivariate density of 2 variables for two
  competing forecasts of those variables.  
}
\usage{
mountains(fcasts1, fcasts2, varnames, pts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fcasts1}{ \eqn{gibbs \times 2}{gibbs x 2} set of forecasts from model 1 }
  \item{fcasts2}{ \eqn{gibbs \times 2}{gibbs x 2} set of forecasts from model 2}
  \item{varnames}{ \code{c("name1","name2")} object of the variable names}
  \item{pts}{ \code{c(pt1,pt2)} which are reference points to be plotted.}
  \item{\dots}{ Other graphics parameters.}
}
\details{
  A "mountain plot" provide a \eqn{2 \times 2}{2 x 2} graph of plots that summarize the
  bivariate forecasts for two competing forecasts. This function
  presents four perspectives on the bivariate density or 'hills' for a
  set of forecasts.  Starting from the bottom right plot and working
  counter-clockwise, the first plot is the bivariate density of the two
  competing forecasts.  The next plot is a contour map that provide the
  topography of the densities.  The third and fourth plots are
  projections of densities in each variable.  The first forecast in the
  function is presented in black, the second in red.  The densities are
  estimated from the Gibbs Monte Carlo sample of forecasts using the
  \code{bkde2D} bivariate kernel density estimator with an optimal
  plug-in bandwidth selected using \code{dpill}.
}
\value{
  None.  Produces the mountain plot described above in the current
  graphics device.
}
\references{ Brandt, Patrick T. and John R. Freeman. 2006. "Advances in
  Bayesian  Time Series Modeling and the Study of Politics: Theory
  Testing, Forecasting, and Policy Analysis"
  \emph{Political Analysis} 14(1):1-36.
}
\author{ Patrick T. Brandt}
\note{ This function requires the bivariate kernel smoother in the
  package \code{\link[KernSmooth]{bkde2D}}}

\seealso{ \code{\link[KernSmooth]{bkde2D}} for details
  of the density estimators
  }
  \examples{
\dontrun{
data(IsraelPalestineConflict)

# Fit a BVAR model
fit.BVAR <- szbvar(IsraelPalestineConflict, p=6, z=NULL, lambda0=0.6,
                   lambda1=0.1, lambda3=2, lambda4=0.5, lambda5=0,
                   mu5=0, mu6=0, nu=3, qm=4, prior=0,
                   posterior.fit=FALSE)

# Fit a flat prior / MLE model
fit.FREQ <- szbvar(IsraelPalestineConflict, p=6, z=NULL, lambda0=0.6,
                   lambda1=0.1, lambda3=2, lambda4=0.5, lambda5=0,
                   mu5=0, mu6=0, nu=3, qm=4, prior=2,
                   posterior.fit=FALSE)

# Generate unconditional forecasts for both models
forecast.BVAR <- uc.forecast.var(fit.BVAR, nsteps=2,
                                 burnin=100, gibbs=1000)

forecast.FREQ <- uc.forecast.var(fit.FREQ, nsteps=2,
                                 burnin=100, gibbs=1000)

# Plot the densities for the forecasts in period of the forecast horizon

mountains(forecast.BVAR$forecast[,2,1:2],
  forecast.FREQ$forecast[,2,1:2], varnames=c("I2P","P2I"), pts=c(0,0))
}
}
\keyword{ hplot}
\keyword{ dplot}
\keyword{ smooth}

