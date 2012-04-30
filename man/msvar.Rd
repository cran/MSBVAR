\name{msvar}
\alias{msvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov-switching vector autoregression (MSVAR) estimator}

\description{
  Estimates a Markov-swtiching vector autoregression (MSVAR) model with
  \eqn{h} regimes (states) by maximum likelihood.  The Hamilton
  filtering algorithm is used to estimate the regimes.  The numerical
  optimization to compute the MLE is based on the block-wise algorithm
  of Sims, Waggoner and Zha (2008).
}
\usage{
msvar(Y, p, h, niterblkopt = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ \eqn{T \times m}{T x m} multiple time series object created
    with \code{ts()}}.
  \item{p}{ Lag length, an integer }
  \item{h}{ Number of regimes / states, an integer }
  \item{niterblkopt}{ Number of iterations to allow for the block-wise
  optimization. }
}
\details{ This function computes ML estimates for an MSVAR(p,h) model
  where \eqn{p} is the number of lags and \eqn{h} is the number of
  regimes. The model is estimated using the block-wise algorithm of
  Sims, Waggoner, and Zha (2008).  This ML optimization algorithm splits
  the parameter space of the MSVAR model into separate block components:
  (1) the transition matrix \eqn{Q}, (2) the
  intercepts, (3) the autoregressive coefficients, (4) the error
  covariances.  The algorithm does 4 separate optimizations for each
  \code{niterblkopt} calls.  Each component of the model is optimized
  separately over the \code{niterblkopt} values using separate calls to
  \code{\link{optim}}. Within each \code{\link{optim}} call, Fortran
  code is used to do the work of the filtering algorithm for the regimes
  in the model
}

\value{
  A list of class \code{MSVAR} and the appropriate inputs objects to
  feed the results into subsequent functions like
  \code{\link{gibbs.msbvar}} (though you should use \code{\link{msbvar}}
  and specify a prior!). 

  \item{init.model }{Description of 'comp1'}
  \item{hreg }{Description of 'comp2'}
  \item{Q }{ \eqn{h \times h}{h x h} Markov transition matrix}
  \item{fp }{ \eqn{T \times h}{T x h} Transition probability matrix}
  \item{m }{ Integer, number of equations}
  \item{p }{ Integer, number of lags}
  \item{h }{ Integer, number of regimes}
  \item{llfval }{ Vector of length \code{niterblkopt} }
  \item{DirectBFGSLastSuccess}{ \code{\link{optim}} convergence code
    returned in the last optimization used in the last block-wise
    optimization}
}
\references{
  Hamilton, James. 1989. "A new approach to the economic analysis of
  nonstationary time series and the business cycle." Econmetrica,
  357--384. 
  
  Sims, Christopher A. and Daniel F. Waggoner and Tao
  Zha. 2008. "Methods for inference in large multiple-equation
  Markov-switching models"  Journal of Econometrics 146(2):255--274.
}
\author{ Patrick T. Brandt and Ryan Davis}

\note{ Consult the \code{\link{msbvar}} function for more details on the
  model. This function is only included as a baseline or helper to the
  overall estimation goal of fitting MSBVAR models.
}

\seealso{ \code{\link{msbvar}} for the Bayesian estimator,
  \code{\link{szbvar}} for the Bayesian, non-regime-switching version,
  \code{\link{gibbs.msbvar}} for posterior sampling. 
}


\examples{
\dontrun{
# Simple replication of Hamilton (1989) as in
# Kim and Nelson (1999: 79, 220)

data(HamiltonGDP)
set.seed(1)

m2 <- msvar(HamiltonGDP, p=1, h=2, niterblkopt=20)

# Now plot the filtered probabilities of a recession
# Compare to Kim and Nelson (1999: 79, 220)

fp.rec <- ts(m2$fp[,1], start=tsp(gdp)[1], freq=tsp(gdp)[3])
plot(fp.rec)

}
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ models }% __ONLY ONE__ keyword per line
