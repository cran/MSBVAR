\name{summary.forecast}
\alias{summary.forecast}
%\alias{summary.forecast.VAR}
%\alias{summary.forecast.BVAR}
%\alias{summary.forecast.BSVAR}
\title{ Summary functions for forecasts obtained through VAR / BVAR /
  B-SVAR model objects}

\description{ Prints a summary of the mean and quantile values for the
  forecasts generated through VAR / BVAR / B-SVAR model objects to
  standard output.
}
\usage{ \method{summary}{forecast}(object, probs = c(0.16,0.84), ...)
}
\arguments{
  \item{object}{ Forecast object genrated through fitting a VAR, BVAR,
    or B-SVAR model
    from either \code{forecast.VAR}, \code{forecast.BVAR}, or
    \code{forecast.BSVAR}}
  \item{probs}{vector list of probability range for quantiles.  default:
  c(0.16,0.84) or a 68\% region (approximately one standard deviation on
  each side of the mean)}
\item{...}{ optional arguments (ignored, but included for S3
  consistency)}
}
\details{
  Prints a summary of the mean and quantile values for the forecasts.
}
\value{
  Returns the mean forecast and the specified posterior
  probability interval for the forecasts.
}
\author{ Patrick T. Brandt}
%\references{}
\seealso{\code{\link{summary}}
}
\examples{
\dontrun{
summary(x)
}
}

\keyword{ print}
\keyword{ models}
\keyword{ manip}
