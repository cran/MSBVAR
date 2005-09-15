\name{forecast.var}
\alias{forecast.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate forecasts for fitted VAR objects}
\description{
  Forecasting for VAR objects with structural (endogenous) and exogenous shocks.
}
\usage{
forecast.var(var.obj, nsteps,
             A0 = t(chol(var.obj$mean.S)),
             shocks = matrix(0, nrow = nsteps, ncol = dim(var.obj$ar.coefs)[1]),
             exog.fut = matrix(0, nrow = nsteps, ncol = nrow(exog.coefs)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var.obj}{ Fitted VAR model of the class "VARobject" produced by
    either \code{mlevar} or \code{szbvar}}
  \item{nsteps}{ Number of periods in the forecast horizon}
  \item{A0}{ M X M matrix of the decomposition of the contemporaneous
    endogenous forecast innovations}
  \item{shocks}{ Structural shocks to the VAR model.  These must be
    scaled consistent with the structural identification in A0}
  \item{exog.fut}{ nsteps x number of exogenous variables matrix of the
    future values of exogenous variable shocks.}
}
\details{
  This function computes VAR forecasts for both classical and Bayesian
  VARs.  Users can specify shocks to the system over the forecast
  horizon (both structural and exogenous shocks).  The forecasting model
  is that described by Waggoner and Zha (1999) and can be used to
  construct conditional and unconditional forecasts based on thes
  structural shocks and the contemporaneous decomposition of the
  innovation variance, A0.
}
\value{
  A matrix time series object, ((T + nsteps) x M) of the original series and
    forecasts}
}
\references{Waggoner, Daniel F. and Tao Zha. 1999. "Conditional
  Forecasts in Dynamic Multivariate Models" Review of \emph{Economics and
  Statistics}, 81(4):639-651.}
}
\author{ Patrick T. Brandt}
\note{ The forecasts can be easily plotted by using the \code{window()}
  command to select the appropriate sample-forecast horizon}

\seealso{ \code{\link{reduced.form.var}} and \code{\link{szbvar}} for estimation methods
  that create the elements needed to forecast}

\examples{
data(IsraelPalestineConflict)
Y.sample1 <- window(IsraelPalestineConflict, end=c(2002, 52))
Y.sample2 <- window(IsraelPalestineConflict, start=c(2003,1))

# Fit a BVAR model
fit.bvar <- szbvar(Y.sample1, p=6, lambda0=0.6, lambda1=0.1, lambda3=2,
                   lambda4=0.25, lambda5=0, mu5=0, mu6=0, prior=0)

# Forecast -- this gives back the sample PLUS the forecasts!

forecasts <- forecast.var(fit.bvar, nsteps=nrow(Y.sample2))
forecasts.only <- forecasts[(nrow(Y.sample1)+1):nrow(forecasts),]

# Plot forecasts and actual data
i2p <- ts(cbind(Y.sample2[,1], forecasts.only[,1]), start=c(2003,1), freq=52)
p2i <- ts(cbind(Y.sample2[,2], forecasts.only[,2]), start=c(2003,1), freq=52)

par(mfrow=c(2,1))
plot(i2p, plot.type=c("single"))
plot(p2i, plot.type=c("single"))
}
\keyword{ ts}
\keyword{ regression }
