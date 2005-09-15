\name{dfev}
\alias{dfev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Decompositions of Forecast Error Variance (DFEV) for VAR models }
\description{
  Computes the M dimensional decomposition of forecast error variance
  for a VAR model.  User can specify the decomposition of the
  contemporaneous innovations.
}
\usage{
dfev(var.obj, A0 = t(chol(var.obj$mean.S)), k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var.obj}{ VAR object created from fitting a VAR model using
    either szbvar() or mlevar() }
  \item{A0}{ Decomposition of the contemporaneous error covariance
    matrix.  Default is to use the Cholesky decomposition of the
    residual error covariance matrix. }
  \item{k}{ Number of periods over which to compute the deccomposition. }
}
\details{
  The decomposition of the forecast error variance (DFEV) provides a
  measure of the relationship among forecast errors or impact of shocks
  to a VAR system. It is computed by finding the moving average
  representation (MAR) of the VAR model and then tracing out the path of
  innovations through the system.  For each of the M innovations in a
  VAR, the amount of the variance in these forecast errors or
  innovations is computed and returned in a table.  The table can be
  accessed via the \code{\link{print.dfev}} and \code{\link{summary.dfev}}
  functions.
}
\value{
  Returns a list with 

  \item{errors }{ M x M x K of the percentage of the innovations in
    variable i explained by the other M variables. }
  \item{std.err}{ M x k dimension matrix of the forecast standard errors. }
  \item{names}{Variable names}
}
\references{ }
\author{ Patrick T. Brandt }
\note{ The interpretation of the DFEV depends on the decomposition of
  the contemporaneous residuals.  In the default case of a Cholesky
  decomposition, this means that the ordering of the variables in the
  decomposition determines the effect of each innovation on the
  subsequent DFEVs.  For high correlated series, this will mean that the
  DFEV is not very robust to the ordering.}

\seealso{ See also \code{\link{print.dfev}} and
  \code{\link{summary.dfev}} for a nicely formatted tables
  and an output example}
\examples{
    data(IsraelPalestineConflict)
    varnames <- colnames(IsraelPalestineConflict)
    fitted.BVAR <- szbvar(IsraelPalestineConflict, p=6, z=NULL,
                           lambda0=0.6, lambda1=0.1,
                           lambda3=2, lambda4=0.25, lambda5=0, mu5=0,
                           mu6=0, nu=3, qm=4, prior=0,
                           posterior.fit=FALSE)

    A0 <- t(chol(fitted.BVAR$mean.S))
    dat.dfev <- dfev(fitted.BVAR, A0, 24)

    print.dfev(dat.dfev)
    summary.dfev(dat.dfev)
}
\keyword{ ts }% at least one, from doc/KEYWORDS

