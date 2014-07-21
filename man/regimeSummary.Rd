\name{regimeSummary}
\alias{regimeSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Regime probability summaries and regime duration estimates based
  on MCMC output for MSBVAR models
}
\description{ Provides summary and quantile computations for regime
  probabilities and regime durations based on MSBVAR MCMC output
}
\usage{
regimeSummary(x, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))
}

\arguments{
  \item{x}{ output from \code{gibbs.msbvar}, the MCMC sampler for the
  MSBVAR models }
  \item{quantiles}{ quantiles one wants to compute, as is done in the
  \code{coda} package.  Defaults are as given above.}

}
\details{
  This function is mainly a wrapper for calls to the \code{coda} package
  to summarize the MCMC output for the transition matrix of an MSBVAR
  model estimated from \code{\link{gibbs.msbvar}}.  It adds labels to
  the output so one know which regime is which in the output.  In the
  summary of the transition matrix \eqn{Q's} elements \eqn{q_{ij}}{q(ij)}
  for the transition from regime \eqn{i} to regime \eqn{j}.

  The ergodic regime probabilities are computed for draw \eqn{k} of the
  MSBVAR MCMC sampler as described in Kim and Nelson (1999):

  eta <- solve(rbind(cbind(diag(h-1) - t(Q)[1:(h-1),1:(h-1)], t(Q)[1:(h-1),h]),
               rep(1,h)))%*%matrix(c(rep(0,h-1),1))

  This is the gives \eqn{N2} draws of the ergodic probabilities of being
  in regime \eqn{k}.  These are summarized again using \code{coda}
  functions.

  Finally, the ergodic regime probabilities can be used to estimate
  expected long run regime durations.  For \eqn{\eta_k}{eta(k)} the
  expected regime duration is \eqn{1/(1-\eta_k)}{1/(1-eta(k)}. This
  again is summarized over the \eqn{N2} draws using \code{coda}
  functions.
}
\value{ Invisible list with 3 elements:

  \item{Q.summary }{ Summary and quantiles of the \code{x$Q.sample} draws of the
    transition matrix}
  \item{lrQ }{ Summary and quantiles of the long run or ergodic regime
    probabilities}
  \item{durations }{ Summary and quantiles of the estimated regime
    durations}
}
\references{
  Kim, Chang-Jin and Charles R. Nelson. 1999. State-Space Models with
  Regime Switching: Classical and Gibbs-Sampling Approaches with
  Applications. Cambridge: MIT Press.
}
\author{ Patrick T. Brandt }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gibbs.msbvar}}, \code{\link{plotregimeid}}, \code{\link{msbvar}}
}
\examples{
\dontrun{
regimeSummary(x)
}
}
\keyword{ print}
\keyword{ models}% __ONLY ONE__ keyword per line
\keyword{ manip}