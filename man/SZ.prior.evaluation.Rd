\name{SZ.prior.evaluation}
\alias{SZ.prior.evaluation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sims-Zha Bayesian VAR Prior Specification Search }
\description{
  Estimates posterior and in-sample fit measures for a reduced form vector
  autoregression model with different specifications of the Sims-Zha
  hyperparameters values.  
}
\usage{
SZ.prior.evaluation(Y, p,
                    lambda0, lambda1, lambda3, lambda4, lambda5,
                    mu5, mu6, z = NULL, nu = ncol(Y) + 1, qm,
                    prior = 0, nsteps, y.future)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ T x m matrix of endogenous variables for the VAR}
  \item{p}{ Lag length}
  
  \item{lambda0}{ List of values, e,g, \code{c(0.7, 0.8, 0.9)} in [0,1],
    Overall tightness of the prior (discounting of prior scale). }
  \item{lambda1}{ List of values, e,g, \code{c(0.05, 0.1, 0.2)} in
    [0,1], Standard deviation or tightness of the prior around the AR(1)
    parameters. }
  \item{lambda3}{ List of values, e,g, \code{c(0, 1, 2)} for Lag decay
    (>0, with 1=harmonic) }
  \item{lambda4}{ List of values, e,g, \code{c(0.15, 0.2, 0.5)} for
    Standard deviation or tightness around the intercept [>0] }
  \item{lambda5}{ Single value for the standard deviation or tightness
    around the exogneous variable coefficients [>0]}
  \item{mu5}{ Single value for sum of coefficients prior weight [>=0]}
  \item{mu6}{ Single value for dummy Initial observations or cointegration prior [>=0]}
  \item{z}{ Exogenous variables }
  \item{nu}{ Prior degrees of freedom = m+1}
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
  \item{prior}{ One of three values: 0 = Normal-Wishart prior, 1 =
    Normal-flat prior, 2 = flat-flat prior (i.e., akin to MLE)}
  \item{nsteps}{ Number of periods in the forecast horizon }
  \item{y.future}{ Future values of the series, nsteps x m for computing
    the root mean squared error and mean absolute error for the fit}
}
\details{
  This function fits a series of BVAR models for the combinations of
  \code{lambda0}, \code{lambda1}, \code{lambda3}, and \code{lambda4}
  provided.  For each possible value of these parameters specified, a
  Sims-Zha prior BVAR model is fit, posterior fit measures are computed, and
  forecasts are generated over \code{nsteps}.  These \code{nstep}
  forecasts are then compared to a new set of data in \code{y.future}
  and root mean sqaured error and mean absolute error measures are
  computed. 
}
\value{
  A matrix of the results with columns corresponding to the values of 
  "lambda0", "lambda1", "lambda3", "lambda4", "lambda5", "mu5",
  "mu6", "RMSE", "MAE", "MargLLF","MargPosterior".
}
\references{ Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
  Time Series Modeling and the Study of Politics: Theory Testing,
  Forecasting, and Policy Analysis"
   \emph{Political Analysis} 14(1):1-36.}
\author{ Patrick T. Brandt }
\note{  The matrix of the results can be usefully plotted using the
  \code{lattice} package.  See the example below.
}

\seealso{ \code{\link{szbvar}}}
\examples{
Y <- EuStockMarkets
results <- SZ.prior.evaluation(window(Y, start=c(1998, 1),
                                       end=c(1998,149)),
                               p=3,
                               lambda0=c(1,0.9),
                               lambda1=c(0.1,0.2),
                               lambda3=c(0,1),
                               lambda4=c(0.1,0.25),
                               lambda5=0,
                               mu5=4,
                               mu6=4, z=NULL,
                               nu=ncol(Y)+1, qm=4,
                               prior=0,
                               nstep=20, 
                               y.future=window(Y, start=c(1998,150)))

# Now plot the RMSE and marginal posterior of the data for each of the
# 6 period forecasts as a function of the prior parameters.  This can
# easily be done using a lattice graphic.

library(lattice)

attach(as.data.frame(results))
dev.new()
xyplot(RMSE ~ lambda0 | lambda1 + lambda3)
dev.new()
xyplot(logMDD ~ lambda0 | lambda1 + lambda3)
dev.new()
xyplot(LLF ~ lambda0 | lambda1 + lambda3)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate}% __ONLY ONE__ keyword per line
