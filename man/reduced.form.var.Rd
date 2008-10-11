\name{reduced.form.var}
\alias{reduced.form.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of a reduced form VAR model}
\description{
  Estimates a reduced form VAR using equation-by-equation seemingly
  unrelated regression (SUR).
}
\usage{
reduced.form.var(Y, p, z=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ \eqn{T \times m}{T x m} multiple time
    series object created with \code{ts()}.}
  \item{p}{ Lag length}
  \item{z}{ \eqn{T \times k}{T x k} exogenous variables in a matrix of
    \eqn{T} rows.  Can be \code{NULL} if there are none.}
}
\details{
  This is a frequentist VAR estimator.  This is a workhorse function ---
  you will want to use other functions such as \code{\link{irf}},
  \code{\link{mc.irf}} or \code{\link{dfev}} to report and interpret the
  results of this object.
}
\value{
  List of class "VAR" with elements,
  \item{intercept}{ Row vector of the \eqn{m} intercepts.}
  \item{ar.coefs}{ \eqn{m \times m \times p}{m x m x p} array of the AR
    coefficients.  The first \eqn{m \times m}{m x m} array is for lag 1,
    the p'th array for lag p.}
  \item{Bhat }{ \eqn{(mp + k + 1) \times m}{(mp + k + 1) x m} matrix of
    the coefficients, where the columns correspond to the variables in
    the VAR.  Intercepts follow the AR coefficients, etc.}
  \item{exog.coefs}{ \eqn{k \times m}{k x m} matrix of exogenous coefficients,
    or \code{NA} if \code{z=NULL}}
  \item{vcv}{ \eqn{m \times m}{m x m} matrix of the maximum likelihood estimate of the
    residual covariance}
  \item{mean.S}{ \eqn{m \times m}{m x m}matrix of the posterior residual covariance.}
  \item{hstar}{ \eqn{mp \times mp}{mp x mp} right hand side
    variables crossproduct.}
  \item{X}{ Right hand side variables for the estimation of BVAR}
  \item{Y}{ Left hand side variables for the estimation of BVAR}
  \item{y}{ Input data (Y)}
}
\references{ Sims, C.A. 1980.  "Macroeconomics and Reality"
  \emph{Econometrica} 48(1): 1-48.
  }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ See also \code{\link{szbvar}} for BVAR models with the Sims-Zha
  prior and \code{\link{szbsvar}} for Bayesian SVAR models with the
  Sims-Zha prior.}

\examples{
data(IsraelPalestineConflict)
rf.var <- reduced.form.var(IsraelPalestineConflict, p=6)
plot(irf(rf.var, nsteps=12))

}
\keyword{ ts}
\keyword{ models}

