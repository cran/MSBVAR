\name{rdirichlet}
\alias{rdirichlet}
\alias{ddirichlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random draws from and density for Dirichlet distribution }
\description{
  Generate draws from a random Dirichlet distribution or compute its density.
}
\usage{
rdirichlet(n, alpha)
ddirichlet(x, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Number of draws }
  \item{alpha}{ Scale matrix, \eqn{h \times h}{h x h}}
  \item{x}{ value to compute density}
}
\details{
 Draws \eqn{n} values for an \eqn{h \times h}{h x h} Dirichlet random
 variable or computes the density.
}
\value{
  \item{x}{ An \eqn{n \times h}{n x h} matrix of the draws}
}
\references{ MCMCpack and gregmisc }

\note{ Based on code from Kevin Quinn in the MCMCpack package and the
  gregmisc package. 
}

\examples{
rdirichlet(2, matrix(rep(1, 4), 2, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
