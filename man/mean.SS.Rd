\name{mean.SS}
\alias{mean.SS}
\alias{sum.SS}
%\alias{var.SS}
\alias{plot.SS}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary measures and plots for MS-B(S)VAR state-spaces}
\description{
  Provides a summary and plotting methods for the \code{SS} class
  objects produced from sampling the posterior of an MSBVAR model.
  These functions provide the mean regime probabilities and a plotting
  method for them. 
}
\usage{
\method{mean}{SS}(x, ...)
\method{sum}{SS}(x, ...)
\method{plot}{SS}(x, ylab="State Probabilities", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{SS class object produced by sampling the posterior of a
    Markov-switching BVAR model in MSBVAR.  These are produced by
    \code{\link{gibbs.msbvar}}. }
  \item{ylab}{y-axis label for the regime plot. }
  \item{\dots}{Other argument or graphics parameters for \code{plot}. }
}
\details{
  The first two provide the sum and mean of the number of
  time periods in each state of Markov-process.  The last produces a
  time series plot of the regime or state probabilities. These are
  computed from the Markov Chain Monte Carlo sample computed from 
  \code{\link{gibbs.msbvar}}
}
\value{
 Mean and sum are \eqn{T \times h}{T x h} matrices for the first
 two summary functions.  The plot function generates a plot in the
 current device.  These are the posterior probability measures of the
 Markov process regimes across \eqn{T} periods.
}
%\references{
%  }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ \code{\link{gibbs.msbvar}}, \code{\link{msbvar}}}
% \examples{
% \dontrun{
% }
% }
\keyword{ models}
\keyword{ hplot}

