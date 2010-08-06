\name{summary}
\alias{summary}
%\alias{summary.VAR}
%\alias{summary.BVAR}
%\alias{summary.BSVAR}
\title{ Summary functions for VAR / BVAR / B-SVAR model objects}
\description{ Prints a summary of the coefficient matrices for various
  VAR / BVAR / B-SVAR model objects to standard output.
}
\usage{ summary(object, ...)
%	summary.VAR(object, ...)
%        summary.BVAR(object, ...)
%        summary.BSVAR(object, ...)
}
\arguments{
  \item{object}{ Fitted VAR, BVAR, or B-SVAR model from either
    \code{reduced.form.var}, \code{szbvar}, or \code{szbsvar}
  }
  \item{...}{other arguments}
  }
\value{ None.
}
\details{
  Prints (posterior) coefficient matrices for each lag and error covariance
  summaries as appropriate.
}
%\references{}
\seealso{\code{\link{summary}}
}
\examples{
\dontrun{
summary(x)
}
}
\author{ Patrick T. Brandt}
\keyword{ print}
\keyword{ models}
\keyword{ manip}
