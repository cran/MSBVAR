\name{print.dfev}
\alias{print.dfev}
\alias{summary.dfev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Printing DFEV tables }
\description{
  Prints decomposition of forecast error variance tables 
}
\usage{
\method{print}{dfev}(x, latex = F, file = NULL, ...)

\method{summary}{dfev}(object, latex = F, file = NULL, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ DFEV object created by \code{\link{dfev}}}
  \item{object}{DFEV object created by \code{\link{dfev}}}
  \item{latex}{ Logical.  T = format results in LaTeX tables,  default
    is F, text output }
  \item{file}{ File for the results.  If \code{NULL}, prints to standard
    output device.}
  \item{...}{Other print and summary arguments}
}
\details{
  Prints DFEV results in a table using \code{xtable} functions.
}
\value{
  None.  
}
%\references{ }
\author{ Patrick T. Brandt}
\seealso{ See \code{\link{dfev}} for an example.}
%\examples{}
\keyword{ print}
\keyword{ models}
\keyword{ manip}
