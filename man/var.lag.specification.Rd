\name{var.lag.specification}
\alias{var.lag.specification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Automated VAR lag specification testing
}
\description{
  Estimates a series of test statistics and measures for VAR lag length selection.
}
\usage{
var.lag.specification(y, lagmax = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ T x m multiple time series}
  \item{lagmax}{ Maximum lag order to be evaluated. Function will return
    lag length tests for all lag orders less than lagmax.}
}
\details{
  Estimates a series of VAR models for 1 to lagmax and returns a
  sequence of \eqn{\chi^2}{Chi2} tests, AIC, BIC and Schwarz criterion values for
  each lag length.
}
\value{
  Results are printed to standard output (screen or file).  In addition,
  a list of two matrices is returned:

  \item{ldets }{ Lag length, log-determinants, \eqn{\chi^2}{Chi2} tests, and
    p-values for each lag length, compared to the null of the next
    shorter lag length}
  \item{results }{Lag length, AIC, BIC, and Scwarz criteria for each lag length}

}
\references{ }
\author{ Patrick T. Brandt }
\note{ Sizes of p-values are uncorrected for multiple testing.  Use cautiously. }


\seealso{ See Also \code{\link{reduced.form.var}} for frequentist VAR estimation
  and \code{\link{szbvar}} for Bayesian VAR estimation.}
\examples{
    data(IsraelPalestineConflict)
    var.lag.specification(IsraelPalestineConflict, lagmax=12)
}
\keyword{ ts}
\keyword{ htest}
  
