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
msbvar(Y, z=NULL, p, h, lambda0, lambda1, lambda3,
       lambda4, lambda5, mu5, mu6, qm,
       alpha.prior=100*diag(h) + matrix(2, h, h),
       prior=0, max.iter=40, initialize.opt=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{Y}{ A \eqn{T \times m}{T x m} multiple time series object created
    with \code{ts()}.}
  \item{z}{ NOT IMPLEMENTED AT PRESENT: THIS SHOULD BE
    A \eqn{T \times k}{T x k} matrix of exogenous variables.  Can be \code{z
      = NULL} if there are none (the default). }
  \item{p}{ Lag length, aninteger} 
  \item{h}{ Number of regimes / states, an integer}
  \item{lambda0}{ Value in \eqn{[0,1]}, Overall tightness of the prior (discounting
    of prior scale).}
  \item{lambda1}{ Value in \eqn{[0,1]}, Standard deviation or
    tightness of the prior around the AR(1) parameters. }
  \item{lambda3}{Lag decay (\eqn{>0}, with 1=harmonic) }
  \item{lambda4}{ Standard
    deviation or tightness around the intercept \eqn{>0} }
  \item{lambda5}{Standard deviation or tightness around the exogneous variable
    coefficients \eqn{>0}}
  \item{mu5}{ Sum of coefficients prior weight
    \eqn{\ge0}.  Larger values imply difference stationarity.}
  \item{mu6}{Dummy initial observations or drift prior \eqn{\ge0}.  Larger values
    allow for common trends.}
  \item{qm}{ Frequency of the data for lag
    decay equivalence.  Default is 4, and a value of 12 will match the lag
    decay of monthly to quarterly data.  Other values have the same effect
    as "4"}
  \item{alpha.prior}{ Prior for the Dirichlet process for the MS
    process.  Default is \code{100 * diag(h) + matrix(2, h, h)}, but the
    model will be sensitive to this.}
  \item{prior}{ One of three values: 0= Normal-Wishart prior, 1 =
    Normal-flat prior, 2 = flat-flat prior
    (i.e., akin to MLE).  The conjugate prior is the first one, which is the
    default.}
  \item{max.iter}{ Maximum number of iterations for the block
    EM algorithm used to fit an initial guess of the model posterior.
    Default value is 40 iterations.  Larger problems will need more
    iterations. }
  \item{initialize.opt }{ Initial values for the block
    optimization algorithm.  If \code{default=NULL}
    \code{\link{initialize.msbvar}} is called to provide values.  User can
    specify values as long as they conform to the structure produced by
    \code{\link{initialize.msbvar}}. }
}

\details{ This function
  estimates the posterior mode of a reduced form Bayesian Markov-switching
  VAR model.  The MSBVAR mode is estimated using block EM algorithm where
  the blocks are 1) the BVAR regression coefficients for each regime
  (separating optimands for intercepts, AR coefficients, and error
  covariances) and 2) the transition matrix.  Starting values are randomly
  drawn, so a random number seed should be set prior to calling the
  function in order to make the results replicable.

  The steps of the blockwise optimization follow the suggestions of Sims,
  Waggoner, and Zha (2008).  The joint optimization problem is
  partitioned into the following separate blocks.  For each block, a
  separate call to \code{\link{optim}} is made, holding all of the other
  blocks constant:  
  \enumerate{
    \item \emph{Maximize over the intercepts}
    \item \emph{Maximize over the AR(p) coefficients}
    \item \emph{Maximize over the error covariances \eqn{\Sigma}{Sigma}}
    \item \emph{Maximize over the transition matrix \eqn{Q}}
  }
 
 These four blocks are iterated a total of \code{max.iter} times.
 Internal to each block, the state-space filtering algorithm for the
 regime classifications is computing using compiled Fortran code for
 speed.   Despite the use of compiled code, this algorithm can take
 several minutes to compute.

 The user should try multiple starting values and number of iterations
 to ensure convergence.  The algorithm will improve with each step of
 the optimization, although sometimes this can be very incremental
 improvement. 

 The results for posterior sampling via \code{\link{gibbs.msbvar}} will
 be sensitive to the choice of \code{alpha.prior}.  This is the prior
 for the independent Dirichlet process for the MS process.  Note that
 the prior is roughly proportionate to the number of time periods spent
 in each regime, since the estimated MS probabilities map to the
 duration of the regime via \eqn{1/(1-p)} where \eqn{p} is the
 probability of staying in the regime.
 
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

  Kim, Chang-Jin and Charles R. Nelson. 1999. State-Space Models with
  Regime Switching: Classical and Gibbs-Sampling Approaches with
  Applications. Cambridge: MIT Press.
  
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
  package vignette to see how this function is used to setup an MSBVAR
  model. An example is currently in \code{\link{gibbs.msbvar}}.}

\seealso{ \code{\link{gibbs.msbvar}} for the MCMC sampler after using
  this function, \code{\link{szbvar}} for a non-switching, Bayesian VAR
  and more details.}
\examples{
\dontrun{
# Simple replication of Hamilton (1989) as in
# Kim and Nelson (1999: 79, 220)

data(HamiltonGDP)
set.seed(214)

m2 <- msbvar(HamiltonGDP, p=1, h=2,
             lambda0=0.8, lambda1=0.15, lambda3=1, lambda4=0.25,
             lambda5=1, mu5=0, mu6=0, qm=12,
             alpha.prior=c(100, 30)*diag(2) +
             matrix(12, 2, 2), prior=0, max.iter=30,
             initialize.opt=NULL)

# Now plot the filtered probabilities of a recession
# Compare to Kim and Nelson (1999: 79, 220)

fp.rec <- ts(m2$fp[,1], start=tsp(HamiltonGDP)[1],
             freq=tsp(HamiltonGDP)[3])
plot(fp.rec)

}
 }
\keyword{ ts }
\keyword{ models }
