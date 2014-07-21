\name{plot.ms.irf}
\alias{plot.ms.irf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Color plot of MSBVAR impulse response functions}

\description{
  Provides an overplotted, color-coded version of the MSBVAR IRFs plot.
  This is an experimental function using color rather than the separate
  plots produced in \code{plot.mc.irf}
}
\usage{
\method{plot}{ms.irf}(x, method = "Sims-Zha2", component = 1,
            probs = c(0.16, 0.84), varnames = attr(x, "eqnames"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Output of the \code{mc.irf} function for an MSBVAR model via
  \code{\link{gibbs.msbvar}}}
  \item{method}{ Method to be used for the error band construction.
    Default method is to use the eigendecomposition method proposed by
    Sims and Zha.  Defined methods are "Percentile" (error bands are
    based on percentiles specified in \code{probs}), "Normal
    Approximation" (Gaussian approximation for interval of width
    \code{probs}), "Sims-Zha1" (Gaussian approximation with linear
    eigendecomposition), "Sims-Zha2" (Percentiles with
    eigendecomposition for each impulse response function), "Sims-Zha3"
    (Percentiles with eigendecomposition of the full stacked impulse
    responses) }
  \item{component}{ If using one of the eigendecomposition methods, the
    eigenvector \code{component} to be used for the error band
    construction.  Default is the first or largest eigenvector
    component. }
  \item{probs}{ is the width of the error bands.  Default
    is \code{c(0.16, 0.84)} which is a 68\% band that is approximately
    one standard deviation, as suggested by Sims and Zha.}
  \item{varnames}{ List of variable names of length \eqn{m} for labeling
    the impulse responses. Default are the input variable names from the
    relevent estimation method.}
  \item{\dots}{ Other graphics parameters. }
}
\details{ This function plots the output of a Monte Carlo simulation of
  MSBVAR impulse response functions produced by \code{mc.irf}.  The function
 allows the user to choose among a variety of frequentist (normal
 appproximation and percentile) and Bayesian (eigendecomposition)
 methods for constructing error bands around a set of impulse
 responses. Impulses or shocks are in the columns and the rows are the
 responses.  Here the plot colors the responses for each reqime, per the
 R default color pallette for colors \code{1:h}.
}
\value{  The primary reason for this function is to plot impulse responses and
  their error bands.  Secondarily, it returns
  an invisible list of the impulses responses, their error bands, and
  summary measures of the fractions of the variance in the eigenvector
  methods that explain the total variation of each response.

  \item{responses }{Responses and their error bands}
  \item{eigenvector.fractions }{Fraction of the variation in each
    response that is explained by the chosen eigenvectors.  \code{NULL}
    for non-eigenvector methods.}
}
\references{
  Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
  Time Series Modeling and the Study of Politics: Theory Testing,
  Forecasting, and Policy Analysis"
  \emph{Political Analysis} 14(1):1-36.

   Sims, C.A. and Tao Zha. 1999. "Error Bands for Impulse
  Responses." \emph{Econometrica}. 67(5): 1113-1156.
}
\author{ Patrick T. Brandt }
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.mc.irf}} }
\examples{
\dontrun{
data(IsraelPalestineConflict)
m1 <- msbvar(IsraelPalestineConflict, p=1, h=2, lambda0=0.6,
             lambda1=0.1, lambda3=1, lambda4=0.5, lambda5=0,
             mu5=0, mu6=0, qm=12, alpha.prior=matrix(10, 2, 2),
             prior=0, max.iter=20)
m2p <- gibbs.msbvar(m1, N1=1000, N2=10000, permute=FALSE, Sigma.idx=1)

irf2 <- mc.irf(m2p, nsteps=12)
plot.ms.irf(irf2)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}
\keyword{ hplot}% __ONLY ONE__ keyword per line
