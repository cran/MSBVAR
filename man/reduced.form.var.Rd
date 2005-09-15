\name{reduced.form.var}
\alias{reduced.form.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of a reduced form VAR model}
\description{
  Estimates a reduced form VAR using equation-by-equation seemingly
  unrelated regression (SUR).
}
\usage{
reduced.form.var(dat, p, z=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ T x m matrix of multiple time series or a multiple time
    series object.}
  \item{p}{ Lag length}
  \item{z}{ Exogenous variables in a matrix of T rows.  Can be
    \code{NULL} if there are none.}
}
\details{
  Frequentist VAR estimation.  This is a workhorse function --- you will
  want to use other functions to report and interpret the results of
  this object.
}
\value{
  List with elements,
  \item{intercept} {Row vector of the m intercepts}
  \item{ar.coefs} { M x M X P array of the AR coefficients.  The first M
  X M array is for lag 1, the p'th array for lag p.}
  \item{Bhat } { (M*P + dim(z)[2]) x M matrix of the coefficients, where the
    columns correspond to the variables in the VAR}
  \item{exog.coefs} { \code{dim(z)[2]} x M matrix of exogenous coefficients,
    or \code{NA} if \code{z=NULL}}
  \item{vcv} { M x M matrix of the maximum likelihood estimate of the
    residual covariance}
  \item{mean.S} {M x M matrix of the posterior residual covariance.}
  \item{hstar} {(M*P) x M right hand side
    variables crossproduct.}
  \item{X} {Right hand side variables for the estimation of BVAR}
  \item{Y} {Left hand side variables for the estimation of BVAR}
  \item{y} {Input data (dat)}
}
\references{ Sims, C.A. 1980.  "Macroeconomics and Reality"
  \emph{Econometrica} 48(1): 1-48.
  }
\author{ Patrick T. Brandt}
\note{ }

\seealso{ See Also \code{\link{szbvar}} for Bayesian estimation.}
\examples{
data(IsraelPalestineConflict)
rf.var <- reduced.form.var(IsraelPalestineConflict, p=6)
plot(irf.var(rf.var, nsteps=12))
}
}
\keyword{ ts}
\keyword{ models}

