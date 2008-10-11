\name{null.space}
\alias{null.space}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the null space of a matrix }
\description{
  Computes the null space of A for an arbitrary linear system of the
  form Ax = b.
}
\usage{
null.space(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \eqn{m \times n}{m x n} \eqn{A} matrix of a linear system
    \eqn{Ax = b} }
}
\details{
  Computes the null space via singular value decomposition (SVD) of A by
  finding the columns of the SVD of A that correspond to the
  non-singular column vectors that span A.
}
\value{
  Returns an \eqn{m \times q}{m x q} matrix that is the null space,
  where \eqn{q} is the rank of \eqn{A}.
}
%\references{  }
\author{ Patrick T. Brandt }
%\note{ }

\seealso{ \code{\link{svd}} }
%\examples{}
\keyword{ algebra }% at least one, from doc/KEYWORDS
\keyword{ array }% __ONLY ONE__ keyword per line
