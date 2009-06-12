\name{irf.msbsvar}
\Rdversion{1.1}
\alias{irf.msbsvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Monte Carlo Integration / Simulation of Impulse Response Functions for
  an MSBSVAR model.
}
\description{
  Simulates a posterior of impulse response functions (IRF) by
  Monte Carlo integration.  This can handle MSBSVAR models estimated
  with the \code{msbsvar} function.  Simulations of IRFs
  from the Bayesian model utilize the posterior estimates for that model.  
}
\usage{
irf.msbsvar(gibbs, msbsvar, nsteps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gibbs}{Output from the \code{gibbs.msbsvar} function.  This is
    the posterior for the MSBSVAR model
}
  \item{msbsvar}{Fitted \code{msbsvar} posterior mode from a call to
    \code{msbsvar}. 
}
  \item{nsteps}{Integer, number of steps or periods in the impulse
    response horizon.
}
}
\details{
  This is a beta function.  It averages over the posterior for the
  regimes and generates a set of regime-specific impulse responses.
  There is one set of responses produced for each endogenous variable in
  each regime.  The identification of the \eqn{A_0}{A(0)} for the system
  is set up in the \code{msbsvar} function.  See that for details.
}
\value{
  
  MSBSVAR:
  
  An irf.MSBSVAR class object object that is the array of
  impulse response samples for the Monte Carlo samples

  \item{impulse}{A list of \eqn{h} \eqn{draws \times nsteps \times m^2}{
      draws X nsteps X (m*m)} arrays of the impulse responses.  The list
      if ordered by regime.  These are IRFs using the associated B-SVAR
      model in \code{msbsvar} and the posterior
      \eqn{A_0}{A(0)}.}
  }
% \references{
% }
\author{ Patrick T. Brandt }
\note{
Eventually this function will be wrapped into the classing mechanism for
estimation and plotting for IRFs.  This is described in the
\code{\link{mc.irf}} documentation.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{gibbs.msbsvar}}, \code{\link{msbsvar}},
  \code{\link{irf}}, \code{\link{mc.irf}}
}

% \examples{
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ models }
\keyword{ regression }
