\name{posterior.fit}
\alias{posterior.fit}
%\alias{posterior.fit.VAR}
%\alias{posterior.fit.BVAR}
%\alias{posterior.fit.BSVAR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimates the marginal likelihood or log posterior probability for
  BVAR, BSVAR, and MSBVAR models}
\description{
  Computes the marginal log likelihood
  other posterior fit measures for BVAR, BSVAR, and MSBVAR models fit with
  \code{\link{szbvar}}, \code{\link{szbsvar}} and, \code{\link{msbvar}}
  (and their posterior samplers).
}
\usage{
posterior.fit(varobj, A0.posterior.obj=NULL, maxiterbs=500)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ Varies for BVAR, BSVAR, or MSBVAR models.
    For a BVAR model, \code{varobj} = output from a call to
    \code{\link{szbvar}}. For a BSVAR model, \code{varobj} = output from
    a call to \code{\link{szbsvar}}. For MSBVAR models,
    \code{varobj} = output from a call to \code{\link{gibbs.msbvar}}. }
  
  \item{A0.posterior.obj}{MCMC Gibbs object for the B-SVAR model \eqn{A_0}{A(0)}
    from \code{\link{gibbs.A0}} }

  \item{maxiterbs}{ Number of iterations for the bridge sampler for
    computing the marginal likelihood for MSBVAR models}
}
\details{
  Estimates the marginal log likelihood, also known as a log marginal
  data density for the various models.  For the BVAR models, this can
  be computed in closed form.  For the BSVAR models the MCMC data
  augmentation method of Chib (1995) is employed.  For the MSBVAR
  models, the importance sampler, reciprocal importance sampler, and
  bridge sampler methods of Fruwirth-Schnatter (2006) are used.
  Consult these references for details (or look at the source code).

  The computations are done using compiled C++ and Fortran code as of
  version 0.3.0.  See the package source code for details about the
  implementation.
}
\value{
  BVAR:

  A list of the class "posterior.fit.VAR" that includes the following
  elements:
  
  \item{data.marg.llf}{ Log marginal density, the probability of the
    data after integrating out the parameters in the model.}
  \item{data.marg.post}{ Predictive marginal posterior density}
  \item{coefficient.post}{ Contribution to the posterior fit from the pdf
    of the coefficients.}
  
  BSVAR:
  
  A list of the class "posterior.fit.BSVAR" that includes the following
  elements:

  \item{log.prior }{ Log prior probability}
  \item{log.llf }{ \eqn{T \times 1}{T x 1} list of the log probabilities
    for each observation conditional on the parameters.}
  \item{log.posterior.Aplus }{ Log marginal probability of
    \eqn{A_1,\ldots,A_p}{A(1),...,A(p)} conditional on the data and
    \eqn{A_0}{A(0)}}
  \item{log.marginal.data.density }{ Log data density or marginal log
    likelihood, the probability of the data after integrating out the
    parameters in the model.}
  \item{log.marginal.A0k }{ \eqn{m \times 1}{m x 1} list of the log
    probabilities of each column (corresponding to the equations) of
    \eqn{A_0}{A(0)} conditional on the other columns. }

  MSBVAR:

  A list of the class "posterior.fit.MSBVAR" that includes the following
  elements:
  
}
\references{

  Chib, Siddartha. 1995. "Marginal Likelihood from the Gibbs Output."
  \emph{Journal of the American Statistical Association}. 90(432):
  1313--1321. 
  
  Waggoner, Daniel F. and Tao A. Zha.  2003. "A Gibbs sampler for
  structural vector autoregressions" \emph{Journal of Economic Dynamics
  \& Control}. 28:349--366.

  Fruhwirth-Schnatter, Sylvia. 2006. Finite Mixture and Markov Switching
  Models. Springer Series in Statistics New York: Springer.,
  esp. Sections 5.4 and 5.5.

  }
\author{ Patrick T. Brandt and W. Ryan Davis}
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
  \code{\link{gibbs.msbvar}}, and 
  \code{\link{print.posterior.fit}} for a print method.
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
