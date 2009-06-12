\name{msbvar}
\alias{msbvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov-switching Bayesian reduced form vector autoregression
  model setup and posterior mode estimation }
\description{
  Sets up and estimates the posterior mode of a reduced form Markov-switching
  Bayesian vector autoregression model with a Sims-Zha prior.  This is the setup and input
  function for the Gibbs sampler for this model.
}
\usage{
msbvar(y, z = NULL, p, h, lambda0, lambda1, lambda3, lambda4,
       lambda5, mu5, mu6, qm,
       alpha.prior = 100 * diag(h) + matrix(2, h, h),
       prior=0, max.iter = 40)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \eqn{T \times m}{T x m} multiple time
    series object created with \code{ts()}.}
  \item{z}{ NOT IMPLEMENTED AT PRESENT: THIS SHOULD BE A \eqn{T \times
      k}{T x k} matrix of exogenous variables.  Can
    be \code{z = NULL} if there are none (the default). }
  \item{p}{ Lag length, an integer}
  \item{h}{ Number of regimes / states, an integer}
  \item{lambda0}{ \eqn{[0,1]}, Overall tightness of the prior (discounting
    of prior scale).}
  \item{lambda1}{ Lag decay (\eqn{>0}, with 1=harmonic) }
  \item{lambda3}{ Lag decay (\eqn{>0}, with 1=harmonic) }
  \item{lambda4}{ Standard deviation or tightness around the intercept \eqn{>0} }
  \item{lambda5}{ Standard deviation or tightness around the exogneous
    variable coefficients \eqn{>0}}
  \item{mu5}{ Sum of coefficients prior weight \eqn{\ge0}.  Larger
    values imply difference stationarity.}
  \item{mu6}{ Dummy initial observations or drift prior \eqn{\ge0}.
    Larger values allow for common trends.}
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
  \item{alpha.prior}{ Prior for the Dirichlet process for the MS
    process.  Default is \code{100 * diag(h) + matrix(2, h, h)}, but the
    model will be sensitive to this.}
  \item{prior}{ One of three values: 0 = Normal-Wishart prior, 1 =
    Normal-flat prior, 2 = flat-flat prior (i.e., akin to MLE).  The
    conjugate prior is the first one, which is the default.}
  \item{max.iter}{ Maximum number of iterations for the block EM
    algorithm used to fit an initial guess of the model posterior.
    Default value is 40 iterations.  Larger problems will need more
    iterations. }
}
\details{
 This function estimates the posterior mode of a reduced form Bayesian
 Markov-switching VAR model.   The MSBVAR mode is estimated using block
 EM algorithm where the blocks are 1) the MS state-space, 2) the BVAR
 regression step for each regime and 3) the transition matrix.  Starting
 values are randomly drawn, so a random number seed should be set prior
 to calling the function in order to make the results replicable.

 This function should NOT be used for inference, since it only finds the
 posterior mode of the model.  This function is intended to generate
 starting values for the Gibbs sampling of the model.  See
 \code{\link{gibbs.msbvar}} for further details of the Gibbs sampling.
}
\value{
  A list describing the posterior mode of the MSBVAR model and the
 inputs necessary for the subsequent Gibbs sampler.
  \item{init.model }{An object of the class BVAR that describes the
    setup of the model.  See \code{\link{szbvar}} for details.}
  \item{hreg }{A list containing the regime-specific moment matrices,
    VAR coefficients, and error covariances}
  \item{Q }{The \eqn{h \times h}{h x h} Markov transition matrix.}
  \item{fp }{The \eqn{T \times h}{T x h} matrix of the filtered regime
    probabilities.  First column is the first regime, etc.}
  \item{m }{Integer, the number of endogenous variables in the system.}
  \item{p }{Integer, the lag length of the VAR.}
  \item{h }{Integer, the number of regimes in the MS process.}
  \item{alpha.prior }{The \eqn{h \times h}{h x h} matrix for the prior
    for the Dirichlet density for the MS process.}
}

\references{
  Brandt, Patrick T. 2009. "Empirical, Regime-Specific Models
  of International, Inter-group Conflict, and Politics"
  
  Fruhwirth-Schnatter, Sylvia. 2001. "Markov Chain Monte Carlo
  Estimation of Classical and Dynamic
  Switching and Mixture Models". Journal of the American Statistical
  Association. 96(153):194--209.

  Fruhwirth-Schnatter, Sylvia. 2006. Finite Mixture and Markov Switching
  Models. Springer Series in Statistics New York: Springer.
  
  Sims, Christopher A. and Daniel F. Waggoner and Tao
  Zha. 2008. "Methods for inference in large multiple-equation
  Markov-switching models"  Journal of Econometrics 146(2):255--274.
  
  Sims, Christopher A. and Tao A. Zha. 1998. "Bayesian Methods for
  Dynamic Multivariate Models" International Economic Review 39(4):949-968.
  
  Sims, Christopher A. and Tao A. Zha. 2006. "Were There Regime Switches
  in U.S. Monetary Policy?" American Economic Review. 96(1):54--81.
  }

\author{ Patrick T. Brandt}
\note{ Users should consult the reference papers and the (coming)
  package vignette to see how this function is used to setup an MSBVAR model.}

\seealso{ \code{\link{gibbs.msbvar}}, \code{\link{szbvar}}}
% \examples{
% }
\keyword{ ts }
\keyword{ models }
