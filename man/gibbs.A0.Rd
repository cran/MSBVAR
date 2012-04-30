\name{gibbs.A0}
\alias{gibbs.A0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gibbs sampler for posterior of Bayesian structural vector
  autoregression models}
\description{
  Samples from the structural contemporaneous parameter matrix \eqn{A_0}{A(0)}
  of a Bayesian Structural Vector Autoregression (B-SVAR) model.  
}
\usage{
gibbs.A0(varobj, N1, N2, thin=1, normalization="DistanceMLA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ A structural BVAR object created by
    \code{\link{szbsvar}} }
  \item{N1}{ Number of burn-in iterations for the Gibbs sampler (should
    probably be greater than or equal to 1000).}
  \item{N2}{ Number of iterations in the posterior sample. }
  \item{thin}{ Thinning parameter for the Gibbs sampler.}
  \item{normalization}{ Normalization rule as defined in
    \code{\link{normalize.svar}}.  Default is "DistanceMLA"
     as recommended in Waggoner and Zha (2003b).}
}
\details{
  Samples the posterior pdf of an \eqn{A_0}{A(0)} matrix for a Bayesian
  structural VAR using the algorithm described in Waggoner and Zha
  (2003a).   This function is meant to be called after
  \code{\link{szbsvar}}, so one should consult that function
  for further information.  The function draws \code{N2 * thin} draws
  from the sampler and returns the \code{N2} draws that are the
  \code{thin}'th elements of the Gibbs sampler sequence.

  The computations are done using compiled C++ code as of version
  0.3.0.  See the package source code for details about the
  implementation.

}
\value{
  A list of class \code{gibbs.A0} with five elements:
  \item{A0.posterior}{A list of three elements containing the results
    of the \code{N2} \eqn{A_0}{A(0)} draws.  The list contains a vector
    storing all of the draws, the location of the drawn elements in and
    the dimension of \eqn{A_0}{A(0)}. \code{A0.posterior$A0} is a vector
    of length equal to the number of parameters in \eqn{A_0}{A(0)} times N2.
    \code{A0.posterior$struct} is a vector of length equal to the number of
      free parameters in \eqn{A_0}{A(0)} that gives the index positions
      of the elements in \eqn{A_0}{A(0)}. \code{A0.posterior$m} is
      \eqn{m}, an integer, the number of equations in the system.
  }
  \item{W.posterior }{ A list of three elements that describes the
    vectorized \eqn{W}
    matrices that characterize the covariance of the restricted
    parameter space of each column of \eqn{A_0}{A(0)}. 
    \code{W.posterior$W} is a vector of the elements of all the sampled
    \eqn{W} matrices.  \code{W.posterior$W.index} is a cumulative index
    of the elements of
      \eqn{W} that defines how the \eqn{W} matrices for each iteration
      of the sampler are stored in the vector. 
    \code{W.posterior$m} is \eqn{m}, an integer, the number of equations
      in the system.}
  \item{ident }{ \code{ident} matrix from the \code{varobj} of binary
    elements that defined the free and restricted parameters, as
    specified in \code{szbsvar}}
  \item{thin }{ \code{thin} value that was input into the function for
    thinning the Gibbs sampler.}
  \item{N2 }{ \code{N2}, size of the posterior sample.}
}
\references{
  Waggoner, Daniel F. and Tao A. Zha.  2003a. "A Gibbs sampler for
  structural vector autoregressions" \emph{Journal of Economic Dynamics
  \& Control}. 28:349--366.

Waggoner, Daniel F. and Tao A. Zha, 2003b. "Likelihood Preserving
Normalization in Multiple Equation Models" \emph{Journal of Econometrics},
114: 329--347}

\author{ Patrick T. Brandt }
\note{ You must have called / loaded an \code{\link{szbsvar}} object to
  use this Gibbs sampler.} 

\seealso{ \code{\link{szbsvar}} for estimation of the
  posterior moments of the B-SVAR model,

  \code{\link{normalize.svar}} for a discussion of and references on
  \eqn{A_0}{A(0)} normalization. 
  
  \code{\link{posterior.fit}} for computing the
  marginal log likelihood for the model after sampling the posterior,

  and \code{\link{plot}} for a unique density plot of the
  \eqn{A_0}{A(0)} elements.
}

\examples{
# SZ, B-SVAR model for the Levant data
data(BCFdata)
m <- ncol(Y)
ident <- diag(m)
ident[1,] <- 1
ident[2,1] <- 1

# estimate the model's posterior moments
set.seed(123)
model <- szbsvar(Y, p=2, z=z2, lambda0=0.8, lambda1=0.1, lambda3=1,
                 lambda4=0.1, lambda5=0.05, mu5=0, mu6=5,
                 ident, qm=12)

# Set length of burn-in and size of posterior.  These are only an
# example.  Production runs should set these much higher.
N1 <- 1000
N2 <- 1000

A0.posterior.obj <- gibbs.A0(model, N1, N2, thin=1)

# Use coda to look at the posterior.
A0.free <- A02mcmc(A0.posterior.obj)

plot(A0.free)
}
\keyword{ ts }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
