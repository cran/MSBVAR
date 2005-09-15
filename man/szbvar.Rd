\name{szbvar}
\alias{szbvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sims-Zha Bayesian VAR model estimation}
\description{
  Estimation of the Bayesian VAR model for just identified VARs
  described in Sims and Zha (1998)
}
\usage{
szbvar(dat, p, z = NULL, lambda0, lambda1, lambda3, lambda4, lambda5,
        mu5, mu6, nu = ncol(dat)+1, qm = 4, prior = 0, posterior.fit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ T x m Matrix of time series or a multiple time series
    object (e.g., ts(matrix(...))) }
  \item{p}{ Lag length }
  \item{z}{ Exogenous variables }
  \item{lambda0}{ [0,1], Overall tightness of the prior (discounting
  of prior scale). }
  \item{lambda1}{ [0,1], Standard deviation or tightness of the prior around the AR(1) parameters. }
  \item{lambda3}{ Lag decay (>0, with 1=harmonic) }
  \item{lambda4}{ Standard deviation or tightness around the intercept [>0] }
  \item{lambda5}{ Standard deviation or tightness around the exogneous
    variable coefficients [>0]}
  \item{mu5}{ Sum of coefficients prior weight [>=0]}
  \item{mu6}{ Dummy Initial observations or cointegration prior [>=0]}
  \item{nu}{ Prior degrees of freedom = m+1}
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
  \item{prior}{ One of three values: 0 = Normal-Wishart prior, 1 =
    Normal-flat prior, 2 = flat-flat prior (i.e., akin to MLE)}
  \item{posterior.fit}{ logical, F = do not estimate log-posterior fit measures,
    T = estimate log-posterior fit measures.}
}
\details{
  This function estimates the Bayesian VAR (BVAR) model described by
  Sims and Zha (1998).  This BVAR model is based a specification of the
  dynamic simultaneous equation representation of the model.  The prior
  is constructed for the structural parameters.
}
\value{
  Returns a list of multiple elements.  This is a workhorse function to
  get the estimates, so nothing is displayed to the screen.  The
  elements of the list are intended as inputs for the various
  post-estimation functions (e.g., impulse response analyses,
  forecasting, decompositions of forecast error variance, etc.)

  Returns a list of the class "VARobject" with the following elements:

  \item{intercept}{Row vector of the m intercepts}
  \item{ar.coefs}{ M x M X P array of the AR coefficients.  The first M
  X M array is for lag 1, the p'th array for lag p.}
  \item{exog.coefs} {dim(z)[2] x M matrix of the coefficients for any
    exogenous variables}
  \item{Bhat }{ (M*P + dim(z)[2]) x M matrix of the coefficients, where the
    columns correspond to the variables in the VAR}
  \item{vcv}{ M x M matrix of the maximum likelihood estimate of the
    residual covariance}
  \item{vcv.Bh}{Posterior estimate of the parameter covariance that is
    conformable with Bhat. }
  \item{mean.S}{M x M matrix of the posterior residual covariance.}
  \item{St} {M x M matrix of the degrees of freedom times the posterior
    residual covariance.}
  \item{hstar}{(M*P + dim(z)[2]) x M prior precision + right hand side
    variables crossproduct.}
  \item{hstarinv}{(M*P + dim(z)[2]) x M prior covariance crossproduct
    (inv(hstar))}
  \item{H0}{(M*P + dim(z)[2]) x M prior precision for the parameters}
  \item{S0}{M x M prior error covariance}
  \item{residuals }{(T-P) x M matrix of the residuals}
  \item{X}{Right hand side variables for the estimation of BVAR}
  \item{Y}{Left hand side variables for the estimation of BVAR}
  \item{y}{Input data (dat)}
  \item{z}{Exogenous variables matrix}
  \item{p}{Lag length}
  \item{num.exog}{Number of exogenous variables}
  \item{qm}{Value of parameter to match quarterly to monthly lag decay
    (4 or 12)}
  \item{prior.type}{Numeric code for prior type: 0 = Normal-Wishart, 1
    = Normal-Flat, 2 = Flat-Flat (approximate MLE)}
  \item{prior} {List of the prior parameter:
    c(lambda0,lambda1,lambda3,lambda4,lambda5, mu5, mu6, nu)}
  \item{marg.llf}{Value of the in-sample marginal log-likelihood for
    the data.}
  \item{marg.post}{Value of the in-sample marginal log posterior of the data.}
  \item{coef.post}{Value of the marginal log posterior estimate of the
    coefficients.}
}
\references{ Sims, C.A. and Tao Zha. 1998. "Bayesian Methods for Dynamic
  Multivariate Models." \emph{International Economic Review}. 39(4):949-968.

  Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
  Time Series Modeling and the Study of Politics: Theory Testing,
  Forecasting, and Policy Analysis". \emph{Political Analysis}.
}
\author{ Patrick T. Brandt, based on code from Robertson and Tallman and
Sims and Zha.}
\note{ This is a work horse function.  You will probably want to use
  other functions to summarize and report the BVAR results.}


\seealso{ \code{\link{reduced.form.var}}. }
\examples{
  \dontrun{
    data(IsraelPalestineConflict)
    varnames <- colnames(IsraelPalestineConflict)

    fit.BVAR <- szbvar(IsraelPalestineConflict, p=6, z=NULL,
                           lambda0=0.6, lambda1=0.1,
                            lambda3=2, lambda4=0.25, lambda5=0, mu5=0,
                            mu6=0, nu=3, qm=4,
                            prior=0, posterior.fit=FALSE)

    # Draw from the posterior pdf of the impulse responses.
    posterior.impulses <- mc.irf.var(fit.BVAR, nsteps=10, draws=5000)

    # Plot the responses
    plot.mc.irf.var(posterior.impulses, method=c("Sims-Zha2"), component=1,
                     probs=c(0.16,0.84), varnames=varnames)
}

}
\keyword{ ts}
\keyword{ models}
