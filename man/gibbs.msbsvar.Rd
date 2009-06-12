\name{gibbs.msbsvar}
\alias{gibbs.msbsvar}
\alias{gibbs.A0.msbsvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gibbs-Metropolis-Hastings sampler for posterior of a
  Markov-switching Bayesian structural vector autoregression model}
\description{
  Samples the parameters for the structural contemporaneous parameter
  matrix \eqn{A_0}{A(0)}, the structural parameters for the lagged
  VAR(p) coefficients, and the state-space for a Markov-switching
  Bayesian Structural Vector Autoregression (MSBSVAR) model.  
}
\usage{
gibbs.msbsvar(x, N1 = 1000, N2 = 1000,
              tune = matrix(10, x$m, x$h))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of the class MSBSVARsetup created by
    \code{\link{msbsvar}} }
  \item{N1}{ Number of burn-in iterations for the sampler}
  \item{N2}{ Number of iterations in the posterior sample}
  \item{tune}{ \eqn{m \times h}{m x h} matrix of tuning values for the
    Metropolis-Hastings steps used to estimate the b parameters of
    \eqn{A_0(s_t)}{A(0, s(t))}}
}

\details{
  Samples the posterior pdf for the MSBSVAR model described in Sims et
  al (2008).   This sampler is rather involved, so users should consult
  the original paper listed in the references below.  This function is
  means to be called after \code{\link{msbsvar}}, so consult that
  function for further information.  The function returns \code{N2}
  draws of the parameters from the sampler.

  Much of the computational heavy lifting (especially the state-space
  sampling for the \eqn{T \times h}{T x h} regime values) is done in
  compiled C++ code.  Consult the course code for additional details.

  This function is experimental, so use at your own risk.
}

\value{
  A list of the class "MSBSVAR" that contains the N2 draws from the
  posterior.   The list has the following elements:

  \item{b.sample}{ a \code{\link[coda]{mcmc}} class object of \code{N2} draws of
    \code{b}}
  \item{F.sample}{ a \code{\link[coda]{mcmc}} class object of \code{N2} draws of the
    coefficients for the F matrix for the AR dynamics.}
  \item{xi.sample}{ a \code{\link[coda]{mcmc}} class object of \code{N2} draws for
    the \eqn{m \times h}{m x h} values of \eqn{\xi}{xi}}
  \item{Q.sample}{ a \code{\link[coda]{mcmc}} class object for the \code{N2} draws of
    the \eqn{h \times h}{h x h} matrix Q}
  \item{ss.sample}{ a list of length \code{N2} of \code{\link[bit]{bit}} class
    objects that store the samples of the 0-1 matrices for the
    \code{h-1} data for the state-space.  These objects can be accessed
    and summarized using the \code{\link{sum.SS}},
    \code{\link{mean.SS}} and \code{\link{plot.SS}} function class.}
  \item{accept.rate}{ \code{mh} vector of the acceptance rates for the
    Metropolis steps for each equation in each regime}
  \item{A0mode}{ \eqn{m \times m}{m x m} matrix for the initial
    \eqn{A_0(s_t)}{A(0, s(t))} around which the draws of the posterior
    are normalized via Waggoner and Zha's likelihood normalization}
  \item{h}{ integer, number of regimes in the model}
}

\references{
  Sims, C.A., D. Waggoner and T. A. Zha. 2008. "Methods for Inference in
  Large Multiple-Equation Markov-switching Models." \emph{Journal of
  Econometrics}. 
}
\author{ Patrick T. Brandt}
\note{ The storage of this object makes extensive use of the
  \code{\link[bit]{bit}} class / package to reduce memory overhead.  You will
  need to use the summary function in MSBVAR to work with these objects
  (or read the source code). }
\seealso{ \code{\link{msbsvar}}, \code{\link{sum.SS}},
  \code{\link{mean.SS}}, \code{\link{plot.SS}}, \code{\link[coda]{mcmc}}  }

%\examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ ts }
\keyword{ regression }% __ONLY ONE__ keyword per line
