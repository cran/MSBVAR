\name{hc.forecast}
\alias{hc.forecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forecast density estimation of hard condition forecasts for VAR
  models via MCMC }
\description{
  Implements a "hard condition" forecast density estimator for VAR/BVAR/B-SVAR
  models as described in Waggoner and Zha (1999).  A "hard condition"
  forecast is one where the forecast path of one or more variables in a
  VAR is constrained to be an exact value.  The forecast densities are
  estimated as the posterior sample for the VAR model using Markov Chain
  Monte Carlo with data augmentation to account for the uncertainty of
  the forecasts and the parameters.  This function DOES account for
  parameter uncertainty in the MCMC algorithm.
}
\usage{
hc.forecast(varobj, yconst, nsteps, burnin,
            gibbs, exog = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ VAR object produced for an unrestricted VAR or BVAR
    using \code{szbvar} or \code{reduced.form.var}}
  \item{yconst}{ nsteps x m matrix of the constrained forecasts that
    matches the variables in the endogenous variables of the VAR object.
   Unconstrained forecasts should be set to NA or zero.}
  \item{nsteps}{ Number of periods in the forecast horizon }
  \item{burnin}{ Burnin cycles for the MCMC algorithm }
  \item{gibbs}{ Number of cycles of the Gibbs sampler after the
    \code{burnin} that are returned in the output}
  \item{exog}{ num.exog x nsteps matrix of the exogenous variable values
  for the forecast horizon.  If left at the \code{NULL} default, they
  are set to zero.}
}
\details{
  "Hard conditions" are restrictions of the future forecast path of a
  variable in a VAR.  Once a variable has been constrained along the
  forecast path, the paths of the other variables in the VAR forecasts
  must be re-estimated to satisfy the forecast constraint, since the
  constrained variable has a forecast variance of zero (it is assumed
  known).  Thus, an MCMC algorithm must be used to determine the
  posterior of the forecasts and a consistent set of VAR parameter
  estimates that satisfy the forecast constraints.  This function
  accounts for the uncertainty of the VAR parameters by sampling from
  them in the computation of the VAR forecasts.
}
\value{
  A list with two components:

  \item{forecast }{ gibbs x nsteps x m array of the samples of the VAR
    forecasts }
  \item{orig.y }{ T x m time series object of the original endogenous
    variables }
}
\references{
   Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
   Time Series Modeling and the Study of Politics: Theory Testing,
   Forecasting, and Policy Analysis"
   \emph{Political Analysis} 14(1):1-36.    

   Waggoner, Daniel F. and Tao Zha. 1999. "Conditional
   Forecasts in Dynamic Multivariate Models" \emph{Review of Economics and
   Statistics}, 81(4):639-651.}
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ \code{\link{plot.forecast}} for plotting,
  \code{\link{forecast}} for unconditional forecasting of forecast
  means, \code{\link{uc.forecast}} for MCMC estimation of forecast
  densities for unconstrained or unconditional forecasts}
\examples{
\dontrun{
## Uses the example from Brandt and Freeman 2006.  Will not run unless
## you have their data from the Politcal
## Analysis website!
library(MSBVAR)   

# Read the data and set up as a time series
data <- read.dta("levant.weekly.79-03.dta") 
attach(data)

# Set up KEDS data
KEDS.data <- ts(cbind(a2i,a2p,i2a,p2a,i2p,p2i),
                start=c(1979,15),
                freq=52,
                names=c("A2I","A2P","I2A","P2A","I2P","P2I"))

# Select the sample we want to use.
KEDS <- window(KEDS.data, end=c(1988,50))


#############################
# Estimate the BVAR models 
#############################

# Fit a flat prior model
KEDS.BVAR.flat <- szbvar(KEDS, p=6, z=NULL, lambda0=1,
                         lambda1=1, lambda3=1, lambda4=1, lambda5=0,
                         mu5=0, mu6=0, nu=0, qm=4, prior=2,
                         posterior.fit=F)

# Reference prior model -- Normal-IW prior pdf
KEDS.BVAR.informed <- szbvar(KEDS, p=6, z=NULL, lambda0=0.6,
                             lambda1=0.1, lambda3=2, lambda4=0.5,
                             lambda5=0, mu5=0, mu6=0,
                             nu=ncol(KEDS)+1, qm=4, prior=0,
                             posterior.fit=F)

# Set up conditional forecast matrix conditions
nsteps <- 12
a2i.condition <- rep(mean(KEDS[,1]) + sqrt(var(KEDS[,1])) , nsteps)

yhat<-matrix(c(a2i.condition,rep(0, nsteps*5)), ncol=6)

# Set the random number seed so we can replicate the results.
set.seed(11023)

# Conditional forecasts
conditional.forcs.ref <- hc.forecast(KEDS.BVAR.informed, yhat, nsteps,
                            burnin=3000, gibbs=5000, exog=NULL)

conditional.forcs.flat <- hc.forecast(KEDS.BVAR.flat, yhat, nsteps,
                             burnin=3000, gibbs=5000, exog=NULL)

# Unconditional forecasts
unconditional.forcs.ref <-uc.forecast(KEDS.BVAR.informed, nsteps,
                                          burnin=3000, gibbs=5000)

unconditional.forcs.flat <- uc.forecast(KEDS.BVAR.flat, nsteps,
                                            burnin=3000, gibbs=5000)

# Set-up and plot the unconditional and conditional forecasts.  This
# code pulls for the forecasts for I2P and P2I and puts them into the
# appropriate array for the figures we want to generate.
uc.flat <- NULL
hc.flat <- NULL
uc.ref <- NULL
hc.ref <- NULL

uc.flat$forecast <- unconditional.forcs.flat$forecast[,,5:6]
hc.flat$forecast <- conditional.forcs.flat$forecast[,,5:6]
uc.ref$forecast <- unconditional.forcs.ref$forecast[,,5:6]
hc.ref$forecast <- conditional.forcs.ref$forecast[,,5:6]

par(mfrow=c(2,2), omi=c(0.25,0.5,0.25,0.25)) 
plot(uc.flat,hc.flat, probs=c(0.16, 0.84), varnames=c("I2P", "P2I"),
     compare.level=KEDS[nrow(KEDS),5:6], lwd=2)
plot(hc.ref,hc.flat, probs=c(0.16, 0.84), varnames=c("I2P", "P2I"),
     compare.level=KEDS[nrow(KEDS),5:6], lwd=2)

}
}
\keyword{ ts}

