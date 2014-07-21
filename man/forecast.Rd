\name{forecast}
\alias{forecast}
%\alias{forecast.VAR}
%\alias{forecast.BVAR}
%\alias{forecast.BSVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate forecasts for fitted VAR objects}
\description{
  Forecasting for VAR/BVAR/BSVAR/MSBVAR objects with structural (endogenous)
  and exogenous shocks.
}
\usage{
forecast(varobj, nsteps, A0=t(chol(varobj$mean.S)),
         shocks=matrix(0, nrow=nsteps, ncol=dim(varobj$ar.coefs)[1]),
         exog.fut=matrix(0, nrow=nsteps, ncol=nrow(varobj$exog.coefs)),
         N1, N2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ Fitted VAR model of the class VAR, BVAR, BSVAR, or MSBVAR
    produced by \code{\link{reduced.form.var}}, \code{\link{szbvar}},
    \code{\link{szbsvar}} or \code{\link{gibbs.msbvar}}.}
  \item{nsteps}{ Number of periods in the forecast horizon}
  \item{A0}{ \eqn{m \times m}{m x m} matrix of the decomposition of the
    contemporaneous endogenous forecast innovations for BSVAR models.}
  \item{shocks}{ Structural shocks to the VAR, BVAR, or BSVAR models.
    These must be scaled consistent with the structural identification
    in \eqn{A_0}{A(0)}.}
  \item{exog.fut}{ nsteps x number of exogenous variables matrix of the
    future values of exogenous variable shocks.  Only implemented for
    VAR, BVAR, and BSVAR models at present.}
  \item{N1}{ integer, number of burnin draws for the MSBVAR forecasts.}
  \item{N2}{ integer, number of final posterior draws for MSBVAR
    forecasts.}
}
\details{
  VAR / BVAR / BSVAR models:

  This function computes forecasts for the classical and Bayesian
  VAR models that are estimated in the MSBVAR package.  Users can
  specify shocks to the system over the forecast horizon (both
  structural and exogenous shocks) for VAR, BVAR, and
  BSVAR models.  The forecasting model is that described by Waggoner and
  Zha (1999) and can be used to construct unconditional forecasts based
  on the structural shocks and the contemporaneous decomposition of the
  innovation variance, A0.

  MSBVAR:

  Generates a set of \code{N2} draws from the posterior forecast
  density.  Forecasts are constructed using data augmentation, so the
  forecasts account for both forecast and parameter uncertainty. The
  function for the MSBVAR model takes as arguments \code{varobj}, which
  is the posterior parameters from a call to \code{\link{gibbs.msbvar}},
  and \code{N1} and \code{N2} to set the burnin and number of draws from
  the posterior.  The posterior forecasts are based on the mixture over
  the \eqn{h} regimes for the specified model.
}
\value{
  For VAR, BVAR, and BSVAR models:

  A matrix time series object,
    \eqn{((T + nsteps) \times m)}{(T + nsteps) x m} of the original
    series and forecasts.

  For MSBVAR models, a list of 4 elements:

  \item{forecasts}{ \eqn{N2 \times nsteps \times m}{N2 x nsteps x m}
    array of the posterior forecasts.}
  \item{ss.sample}{ bit compressed version of the MS state space. (can
    be summarized with \code{\link{plot.SS}} or \code{\link{mean.SS}}.)}
  \item{k}{ number of forecast steps, nsteps}
  \item{h}{ integer, number of MS regimes used in the forecasts.}
}
\references{Waggoner, Daniel F. and Tao Zha. 1999. "Conditional
  Forecasts in Dynamic Multivariate Models" Review of \emph{Economics and
  Statistics}, 81(4):639-651.
}
\author{ Patrick T. Brandt}
\note{ The forecasts can be plotted using the \code{plot.forecast()}
  command to select the appropriate sample-forecast horizon.}

\seealso{ \code{\link{reduced.form.var}}, \code{\link{szbvar}} and
  \code{\link{szbsvar}} for estimation methods that create the elements
  needed to forecast}

\examples{
data(IsraelPalestineConflict)
Y.sample1 <- window(IsraelPalestineConflict, end=c(2002, 52))
Y.sample2 <- window(IsraelPalestineConflict, start=c(2003,1))

# Fit a BVAR model
fit.bvar <- szbvar(Y.sample1, p=6, lambda0=0.6, lambda1=0.1, lambda3=2,
                   lambda4=0.25, lambda5=0, mu5=0, mu6=0, prior=0)

# Forecast -- this gives back the sample PLUS the forecasts!

forecasts <- forecast(fit.bvar, nsteps=nrow(Y.sample2))
forecasts.only <- forecasts[(nrow(Y.sample1)+1):nrow(forecasts),]

# Plot forecasts and actual data
i2p <- ts(cbind(Y.sample2[,1], forecasts.only[,1]),
          start=c(2003,1), freq=52)

p2i <- ts(cbind(Y.sample2[,2], forecasts.only[,2]),
          start=c(2003,1), freq=52)

par(mfrow=c(2,1))
plot(i2p, plot.type=c("single"))
plot(p2i, plot.type=c("single"))


\dontrun{
# MSBVAR forecasts

# Fit model
m1 <- msbvar(Y.sample1, p=1, h=2, lambda0=0.8, lambda1=0.2,
             lambda3=1, lambda4=0.2, lambda5=0, mu5=0, mu6=0,
             qm=12, prior=0)

# Gibbs sampling
m1id <- gibbs.msbvar(m1, N1=1000, N2=10000, permute=FALSE, Sigma.idx=1)

# Forecast density estimation
msforc <- forecast(m1id, nsteps=nrow(Y.sample2), N1=1000, N2=10000)

# Summarize forecasts
apply(msforc$forecasts, c(2,3), mean)

}
}
\keyword{ ts}
\keyword{ regression }
