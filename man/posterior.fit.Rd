\name{posterior.fit}
\alias{posterior.fit}
%\alias{posterior.fit.VAR}
%\alias{posterior.fit.BVAR}
%\alias{posterior.fit.BSVAR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimates the marginal likelihood and posterior probability for
  VAR, BVAR, and BSVAR models}
\description{
  Computes the marginal log likelihood
  other posterior fit measures for VAR, BVAR, and BSVAR models fit with
  \code{\link{szbsvar}} and \code{\link{szbvar}}.
}
\usage{
posterior.fit(varobj, A0.posterior.obj=NULL)
%posterior.fit.VAR(varobj)
%posterior.fit.BVAR(varobj)
%posterior.fit.BSVAR(varobj, A0.posterior.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ BVAR or BSVAR object from fitting a model with
    \code{\link{szbsvar}} or \code{\link{szbvar}} }
  
  \item{A0.posterior.obj}{MCMC Gibbs object for the B-SVAR model \eqn{A_0}{A(0)}
    from \code{\link{gibbs.A0}} }
}
\details{
  Estimates the marginal log likelihood, log prior, log posterior for
  the \eqn{A_0} and \eqn{A_1,\ldots,A_p} parameters of the model, and
  the log data density for the sample (after integrating out the model
  parameters).  The approach used is that of Chib (1995).

  The computations are done using compiled C++ code as of version
  0.3.0.  See the package source code for details about the
  implementation.

}
\value{
  BVAR:

  A list of the class "posterior.fit.VAR" that includes the following
  elements:
  
  \item{data.marg.llf}{Log marginal density, the probability of the
    data after integrating out the parameters in the model.}
  \item{data.marg.post}{Predictive marginal posterior density}
  \item{coefficient.post}{Contribution to the posterior fit from the pdf
    of the coefficients.}
  
  BSVAR:
  
  A list of the class "posterior.fit.BSVAR" that includes the following
  elements:

  \item{log.prior }{ Log prior probability}
  \item{log.llf }{ \eqn{T \times 1}{T x 1} list of the log probabilities
    for each observation conditional on the parameters.}
  \item{log.posterior.Aplus }{ Log marginal probability of
    \eqn{A_1,\ldots,A_p}{A(1),...,A(p)} conditional on the data and \eqn{A_0}{A(0)}}
  \item{log.marginal.data.density }{ Log data density or marginal log
    likelihood, the probability of the data after integrating out the
    parameters in the model.}
  \item{log.marginal.A0k }{ \eqn{m \times 1}{m x 1} list of the log probabilities of each
    column (corresponding to the equations) of \eqn{A_0}{A(0)} conditional on
    the other columns. }
}
\references{

  Chib, Siddartha. 1995. "Marginal Likelihood from the Gibbs Output."
  \emph{Journal of the American Statistical Association}. 90(432):
  1313--1321. 
  
  Waggoner, Daniel F. and Tao A. Zha.  2003. "A Gibbs sampler for
  structural vector autoregressions" \emph{Journal of Economic Dynamics
  \& Control}. 28:349--366.

  }
\author{ Patrick T. Brandt }
\note{ The log Bayes factor for two model can be computed using the
  log.marginal.data.density:

  log BF = log.marginal.data.density.1 - log.marginal.data.density.2

  Note that at present, the scale factors for the BVAR and B-SVAR models
  are different (one used the concentrated likelihood, the other does
  NOT).  Thus, one cannot compare fit measures across the two
  functions.  To compare a recursive B-SVAR to a non-recursive B-SVAR
  model, one should estimate the recursive model with \code{szbsvar}
  using the appropriate \code{ident} matrix and then call
  \code{posterior.fit} on the two B-SVAR models!
  
}

\seealso{ 
  \code{\link{szbvar}},
  \code{\link{szbsvar}},
  \code{\link{gibbs.A0}},
  \code{\link{mc.irf}},
  \code{\link{print.posterior.fit}}
}

\examples{
\dontrun{
varobj <- szbsvar(Y, p, z = NULL, lambda0, lambda1, lambda3, lambda4,
                  lambda5, mu5, mu6, ident, qm = 4)
A0.posterior <- gibbs.A0(varobj, N1, N2)
fit <- posterior.fit(varobj, A0.posterior)
print(fit)
}
}

\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ ts }% __ONLY ONE__ keyword per line
