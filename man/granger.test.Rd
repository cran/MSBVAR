\name{granger.test}
\alias{granger.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bivariate Granger causality testing}
\description{
  Bivariate Granger causality testing for multiple time series.
}
\usage{
granger.test(y, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ T x m time series or matrix.}
  \item{p}{ Lag length to be used for computing the test }
}
\details{
  Estimates all possible bivariate Granger causality tests for m
  variables.  Bivariate Granger causality tests for two variables X and
  Y evaluate whether the past values of X are useful for predicting Y
  once Y's history has been modeled.  The null hypothesis is that the
  past p values of X do not help in predicting the value of Y.

  The test is implemented by regressing Y on p past values of Y and p
  past values of X.  An F-test is then used to determine whether the
  coefficients of the past values of X are jointly zero.

  This produces a matrix with m*(m-1) rows that are all of
  the possible bivariate Granger causal relations.  The results include
  F-statistics and p-values for each test.  Tests are estimated using
  single equation OLS models.  

}
\value{
A matrix with 2 columns.  Column 1 are the F-statistic values.  Column 2
are the p-values for the F-tests.  Row labels specifying the Granger
causality relationship tested will be included if variables in the input
time series \code{y} include variable or \code{dimnames}.
}
\references{
  Granger, C.W.J. 1969. "Investigating Causal Relations by Econometric
  Models and Cross-Spectral Methods" Econometrica 37:424-438.
  
  Sims, C.A. 1972. "Money, Income, and Causality" American
  Economic Review. 62:540-552.
}  
\author{ Patrick T. Brandt }
\note{These are bivariate tests -- not block exogeneity tests for a
  fitted VAR model.  Note also that these tests are highly sensitive to
  lag length (p) and the presence of unit roots.  Results in the matrix
  include row labels for nice printing with \code{xtable()}}

\seealso{ \code{\link{reduced.form.var}} for frequentist VAR estimation,
  \code{\link{szbvar}} for Bayesian VAR estimation with Sims-Zha prior,
  \code{\link{var.lag.specification}} for VAR lag length testing.
  }
\examples{
data(IsraelPalestineConflict)
granger.test(IsraelPalestineConflict, p=6)

}
\keyword{ ts }% at least one, from doc/KEYWORDS
\keyword{ htest}
