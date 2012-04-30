\name{print.posterior.fit}
\alias{print.posterior.fit}
%\alias{print.posterior.fit.VAR}
%\alias{print.posterior.fit.BVAR}
%\alias{print.posterior.fit.BSVAR}
%\alias{print.posterior.fit.MSBVAR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print method for posterior fit measures}
\description{
  Prints objects of the classes "posterior.fit.VAR",
  "posterior.fit.BVAR", and "posterior.fit.BSVAR".
}
\usage{
\method{print}{posterior.fit}(x, ...)
%print.posterior.fit.VAR(x, ...)
%print.posterior.fit.BVAR(x, ...)
%print.posterior.fit.BSVAR(x, ...)
%print.posterior.fit.MSBVAR(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object produced by \code{\link{posterior.fit}}},
  \item{...}{ other print options}
}
\details{
  Called for its side effect --- printing the output of
  \code{\link{posterior.fit}}
  
}
\value{
None
}
%\references{ }
\author{ Patrick T. Brandt}
%\note{ }

\seealso{
  \code{\link{szbvar}},
  \code{\link{szbsvar}},
  \code{\link{gibbs.A0}},
  \code{\link{gibbs.msbvar}},
  \code{\link{mc.irf}},
  \code{\link{posterior.fit}}
}
\examples{
\dontrun{
varobj <- szbsvar(Y, p, z = NULL, lambda0, lambda1, lambda3, lambda4,
                  lambda5, mu5, mu6, ident, qm = 4)
A0.posterior <- gibbs.A0(varobj, N1, N2)
fit <- posterior.fit(varobj, A0.posterior)
print(fit)
}

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ ts }% __ONLY ONE__ keyword per line
