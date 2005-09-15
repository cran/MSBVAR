\name{mc.irf.var}
\alias{mc.irf.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Monte Carlo Integration / Simulation of Impulse Response
  Functions }
\description{
  Simulates a posterior sample of impulse response functions (IRF) by
  Monte carlo integration.  This can handle both Bayesian and
  frequentist VARs estimated with the szbvar() and mlevar() functions.
  The decomposition of the contemporaneous innovations is handled by a
  Cholesky decomposition of the error covariance matrix in each VAR
  object.  Simulations of IRFs from the Bayesian model utilize the
  posterior estimates for that model.
}
\usage{
mc.irf.var(varobj, nsteps, draws)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varobj}{ VAR objects for a fitted VAR model from either
    \code{szbvar()} or \code{reduced.form.var()}}
  \item{nsteps}{ Number of periods over which to compute the impulse responses}
  \item{draws}{ Number of draws for the simulation of the posterior
    distribution of the IRFs }
}
\details{
  Draws a set of posterior samples from the VAR coefficients and
  computes impulse responses for each sample.  These samples can then be
  summarized to compute MCMC based estimates of the responses using the
  error band methods described in Sims and Zha (1999).
}
\value{

  An mc.irf.var class object object that is the array of impulse
  response samples for the Monte Carlo samples

  \item{impulse}{\eqn{draws \times nsteps \times m^2}{ draws X nsteps X
      (m*m)} array of the impulse responses }
}
\references{
  Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
   Time Series Modeling and the Study of Politics: Theory Testing, Forecasting, and Policy Analysis"
   \emph{Political Analysis}.
   
  Sims, C.A. and Tao Zha. 1999. "Error Bands for Impulse
  Responses." \emph{Econometrica} 67(5): 1113-1156.

  Hamilton, James. 1994.  Time Series Analysis.  Chapter 11.
  }
\author{ Patrick T. Brandt }
\note{ }


\seealso{ See also as \code{\link{plot.mc.irf.var}} for plotting methods
  and error band construction for the posterior of the impulse response functions}
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
