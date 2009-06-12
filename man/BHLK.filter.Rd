\name{BHLK.filter}
\alias{BHLK.filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Baum-Hamilton-Lindgren-Kim state-space filter}
\description{
  Implements a Baum-Hamilton-Lindgren-Kim state-space filter for a
  multivariate Markov-switching model.
}
\usage{
BHLK.filter(u, Sigma, Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u }{\eqn{T \times m \times h}{T x m x h} dimensional array of
    the regime-specific residuals to be filtered for the MS process.}
  \item{Sigma }{\eqn{m \times m \times h}{m x m x h} dimensional array
    of the regime-specific error covariances.}
  \item{Q }{The \eqn{h \times h}{h x h} Markov transition matrix.}
}
\details{
  Estimates the Markov-switching regime probabilities for an MSBVAR
  model.  Conditional on the BVAR(p) model, the residuals are used to
  filter the data.  The estimation is done using compiled C++ code for
  speed. 
}
\value{
  \item{x }{ \eqn{T \times h}{T x h} matrix of the filtered regime
    probabilities. }
}
\references{
  Sims, Christopher A. and Daniel F. Waggoner and Tao
  Zha. 2008. "Methods for inference in large multiple-equation
  Markov-switching models"  Journal of Econometrics 146(2):255--274.

  Kim, C.J. and C.R. Nelson. 1999.  State-space models with regime
  switching. Cambridge, Mass:  MIT Press.
  
  Krolzig, Hans-Martin. 1997. Markov-Switching Vector Autoregressions:
  Modeling, Statistical Inference, and Application to Business Cycle
  Analysis. 

}
\author{ Patrick T. Brandt}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{ \code{\link{msbvar}} }
\examples{
# Sets up and estimates a set of filter probabilities
h <- 2
m <- 2
TT <- 100
u <- array(rnorm(TT*m*h), c(TT, m, h))
Sigma <- array(diag(m), c(m, m, h))
Q <- matrix(c(0.99, 0.2, 0.01, 0.8), h, h)
x <- BHLK.filter(u, Sigma, Q)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ models}
