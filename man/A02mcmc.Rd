\name{A02mcmc}
\alias{A02mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Converts A0 objects to coda MCMC objects}

\description{
  Converts A0 objects from \code{gibbs.A0.BSVAR} into
  \code{mcmc} objects for analysis with \code{coda}
}
\usage{
A02mcmc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \eqn{N2 \,\times }{N2 x} number of free parameters in A(0) MCMC Gibbs
    sample object for the B-SVAR model \eqn{A_0}{A(0)}
    from \code{\link{gibbs.A0}}.  This matrix is a column major to row
    major version of A(0) that can be used to diagnose covergence and
    summarize the elements of A(0)}
}
\details{
  Returns an object of the class \code{mcmc}, an N2 x number free
  parameters in A(0) matrix.  This can then be fed into \code{coda} for
  further analysis of the posterior.
}
\value{
  Object with class \code{mcmc}
}

\author{ Patrick T. Brandt }
% \note{ 
% }

\seealso{ 
  \code{\link{gibbs.A0}},\code{\link[coda]{mcmc}}
}



\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ ts }% __ONLY ONE__ keyword per line
