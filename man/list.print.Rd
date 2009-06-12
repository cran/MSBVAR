\name{list.print}
\Rdversion{1.1}
\alias{list.print}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prints a list object for the VAR and BVAR models in MSBVAR
  }

\description{
  Provides a smartly formatted print method for the list objects created
  by MSBVAR objects.  This will provide a table of estimates for the VAR
  and BVAR methods in this package.
}
\usage{
list.print(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Fitted model object from \code{\link{szbvar}} or \code{\link{szbsvar}}
}
}
\details{
  This is a way to view the coefficients from a B(S)-VAr model fit with
  this package.
}
\value{
None.  Results are send to STDOUT.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Patrick T. Brandt and Justin Appleby.
}
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{szbvar}}, \code{\link{szbsvar}}
}
% \examples{
% }

\keyword{ ts }
\keyword{ models }% __ONLY ONE__ keyword per line
