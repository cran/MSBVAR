\name{simulateMSAR}
\alias{simulateMSAR}

\title{ Simulate (univariate) Markov-switching autoregressive (MSAR) data
}
\description{ Simulate (univariate) Markov-switching autoregressive
  (MSAR) data
}
\usage{
simulateMSAR(bigt, Q, theta, st1, y1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bigt}{ Integer, number of observations to generate.
}
  \item{Q}{ \eqn{ h} dimensional transition matrix for the MS
    process. \eqn{h \times h}{h x h} Markov transition matrix whose rows
    sum to 1 with the main weights on the diagonal elements.
}
  \item{theta}{ Matrix of the MSAR coeffients with \eqn{h} rows and
    \eqn{m \times p + 2}{m x p + 2} columns.  The first column
    is the constants, the next \eqn{m \times p + 1}{m x p + 1} columns
    are the autoregressive coefficients (by lag -- so the first \eqn{m
      \times 1}{m x 1} are the AR(1) coefficients, etc.) and the last
    \eqn{m \times 1}{ m x 1} elements are the error variances (remember,
    this is univariate!)
}
  \item{st1}{Starting regime, an integer less than or equal to \eqn{h}
}
  \item{y1}{ Starting value for simulated data in regime \code{st1}
}
}
\details{ This function simulates a univariate MSAR model. The user
  needs to input the transition matrix \eqn{Q} and the autoregression
  coefficients via \eqn{theta}.  The assumption in this model is that
  the error process is Gaussian.
}
\value{
  A list with two elements:
  \item{Y }{ The simulated univariate MSAR time series}
  \item{st }{ A vector of integers identifying the regime of each
    observation in \code{Y}}

}
\references{
  Kim, Chang-Jin and Charles R. Nelson. 1999. State-Space Models with
  Regime Switching: Classical and Gibbs-Sampling Approaches with
  Applications. Cambridge: MIT Press.

}
\author{
  Patrick T. Brandt and Ryan Davis
}
% \note{
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{simulateMSVAR}} for the multivariate version 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example of call here
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ model }% __ONLY ONE__ keyword per line

