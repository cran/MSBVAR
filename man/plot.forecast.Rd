\name{plot.forecast}
\alias{plot.forecast}
%\alias{plot.forecast.VAR}
%\alias{plot.forecast.BVAR}
%\alias{plot.forecast.BSVAR}
\title{Plot function for forecasts}
\description{ Generates simple plots of forecasts obtained from
  forecast.VAR / forecast.BVAR / forecast.B-SVAR}

\usage{ \method{plot}{forecast}(x, ...)
}
\arguments{
  \item{x}{ Plots generated from forecast genrated through fitted VAR,
    BVAR, or B-SVAR model from \code{\link{forecast}}}.
  \item{...}{ Other graphics parameters}
}

\details{
  Generates a plot in the current graphics device for the m time series
  in the respective (B)VAR model.
}


\value{ None.  Generates a plot in the current graphics device.}
%\references{}
\author{ Patrick T. Brandt}
\seealso{\code{\link{summary}}
}
\examples{
\dontrun{
plot(x)
}
}
\keyword{ print}
\keyword{ models}
\keyword{ manip}
