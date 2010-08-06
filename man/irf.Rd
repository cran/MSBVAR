\name{irf}
\alias{irf}
\alias{irf.VAR}
\alias{irf.BVAR}
\alias{irf.BSVAR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Impulse Response Function (IRF) Computation for a VAR}
\description{
  Computes the impulse response function (IRF) or moving average
  representation (MAR) for an m-dimensional set of VAR/BVAR/B-SVAR
  coefficients. 
}
\usage{
irf(varobj, nsteps, A0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{VAR, BVAR, or BSVAR objects for a fitted VAR, BVAR, or BSVAR model
    from \code{szbvar}, \code{szbsvar} or \code{reduced.form.var}}
    \item{nsteps}{ Number or steps, or the horizon over which to compute
    the IRFs (typically 1.5 to 2 times the lag length used in estimation}
  \item{A0}{ Decomposition contemporaneous error covariance of a VAR/BVAR/BSVAR,
    default is a Cholesky decomposition of the error covariance matrix
    for VAR and BVAR models,
    \code{A0 = chol(varobj$mean.S)}, and the inverse of \eqn{A_0}{A(0)}
    for B-SVAR models, \code{A0 = solve(varobj$A0.mode)}}
}
\details{
  This function should rarely be called by the user.  It is a working
  function to compute the IRFs for a VAR model.  Users will typically
  want to used one of the simulation functions that also compute error
  bands for the IRF, such as \code{mc.irf} which calls this function
  and simulates its multivariate posterior distribution.
}
\value{
  A list of the AR coefficients used in computing the IRF and the
  impulse response matrices:
  \item{B}{ \eqn{m \times m \times nstep}{ m x m x nstep} Autoregressive
    coefficient matrices in lag order.  Note that
    all AR coefficient matrices for \eqn{nstep>p}{nstep > p} are zero.}
  \item{mhat }{ \eqn{m \times m \times \times nstep}{m x m x nstep}
    impulse response matrices.  \code{mhat[,,i]} are the
    impulses for the i'th period for the \eqn{m} variables.}
}
\references{ Sims, C.A. and Tao Zha. 1999. "Error Bands for Impulse
  Responses." Econometrica 67(5): 1113-1156.

  Hamilton, James. 1994.  Time Series Analysis.  Chapter 11.
  }
\author{ Patrick T. Brandt}
\note{ The IRF depends on the ordering of the variables and the
  structure of the decomposition in A0.}

\seealso{ See also \code{\link{dfev}} for the related decompositions of
  the forecast error variance, \code{\link{mc.irf}} for Bayesian and
  frequentist computations of IRFs and their variances (which is what
  you probably really want).}
\examples{
data(IsraelPalestineConflict)
rf.var <- reduced.form.var(IsraelPalestineConflict, p=6)
plot(irf(rf.var, nsteps = 12))
}
\keyword{ ts}
\keyword{ models}
