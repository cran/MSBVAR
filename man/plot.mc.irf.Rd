\name{plot.mc.irf}
\alias{plot.mc.irf}
\alias{plot.mc.irf.VAR}
\alias{plot.mc.irf.BVAR}
\alias{plot.mc.irf.BSVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting posteriors of Monte Carlo simulated impulse responses}
\description{
  Provides a plotting method for the \code{mc.irf} Monte Carlo
  sample of impulse responses.  Responses can be plotted with classical
  or Bayesian error bands, as suggested by Sims and Zha (1999).
}
\usage{
plot.mc.irf(x, method=c("Sims-Zha2"), component=1, probs=c(0.16,0.84), 
            varnames = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output of the \code{mc.irf} function}
  \item{method}{ Method to be used for the error band construction.
    Default method is to use the eigendecomposition method proposed by
    Sims and Zha.  Defined methods are "Percentile" (error bands are
    based on percentiles specified in \code{probs}), "Normal
    Approximation" (Gaussian approximation for interval of width
    \code{probs}), "Sims-Zha1" (Gaussian approximation with linear
    eigendecomposition), "Sims-Zha2" (Percentiles with
    eigendecomposition for each impulse response function), "Sims-Zha3"
    (Percentiles with eigendecomposition of the full stacked impulse responses)}
  \item{component}{ If using one of the eigendecomposition methods, the
    eigenvector \code{component} to be used for the error band
    construction.  Default is the first or largest eigenvector component.}
  \item{probs}{ \code{probs} is the width of the error bands.  Default
    is \code{c(0.16, 0.84)} which is a 68\% band that is approximately
    one standard deviation, as suggested by Sims and Zha.}
  \item{varnames}{ List of variable names for labeling the impulse responses.}
  \item{\dots}{ Other graphics parameters}
}
\details{
 This function plots the output of a Monte Carlo simulation of VAR
 impulse response functions produced by \code{mc.irf}.  The function
 allows the user to choose among a variety of frequentist (normal
 appproximation and percentile) and Bayesian (eigendecomposition)
 methods for constructing error bands around a set of impulse 
 responses. Impulses or shocks are in the columns and the rows are the
 responses.
}
\value{

  The primary reason for this function is to plot impulse responses and
  their error bands.  Secondarily, it returns
  an invisible list of the impulses responses, their error bands, and
  summary measures of the fractions of the variance in the eigenvector
  methods that explain the total variation of each response.
  
  \item{responses }{Responses and their error bands}
  \item{eigenvector.fractions }{Fraction of the variation in each
    response that is explained by the chosen eigenvectors.  \code{NULL}
    for non-eigenvector methods.}

}
\references{ Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
  Time Series Modeling and the Study of Politics: Theory Testing,
  Forecasting, and Policy Analysis"
  \emph{Political Analysis} 14(1):1-36.

   Sims, C.A. and Tao Zha. 1999. "Error Bands for Impulse
  Responses." \emph{Econometrica}. 67(5): 1113-1156.

  }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{ See Also \code{\link{mc.irf}} for the computation of Monte
  Carlo samples of impulse responses,
  \code{\link{szbsvar}} for estimation of the posterior
  moments of the B-SVAR model,
  \code{\link{gibbs.A0.BSVAR}} for Gibbs sampling the
  posterior of the \eqn{A_0}{A(0)} for the model, and 
}
\examples{
\dontrun{
data(IsraelPalestineConflict)
fit.BVAR <- szbvar(IsraelPalestineConflict, p=6, z=NULL, lambda0=0.6,
                   lambda1=0.1, lambda3=2, lambda4=0.5, lambda5=0,
                   mu5=0, mu6=0, nu=3, qm=4, prior=0,
                   posterior.fit=FALSE)

posterior.impulses <- mc.irf(fit.BVAR, nsteps=12, draws=1000)
plot(posterior.impulses, method = c("Percentile"))
}
}
\keyword{ models}
\keyword{ hplot}

