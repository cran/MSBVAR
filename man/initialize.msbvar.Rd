\name{initialize.msbvar}
\alias{initialize.msbvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Initializes the mode-finder for a Markov-switching Bayesian VAR model}

\description{ Sets up the initial values for the mode optimization of an
  MSBVAR model with a Sims-Zha prior.  This sets up the
  \code{initialize.opt} argument of the \code{\link{msbvar}} function.
  Users can inputs values outside of the defaults for the \code{Q}
  transition matrix and other arguments with this function.  This
  function also serves as a model for alternative, user-defined initial
  values for the Gibbs sampler.
}
\usage{
initialize.msbvar(y, p, z = NULL, lambda0, lambda1, lambda3, lambda4,
                  lambda5, mu5, mu6, nu, qm, prior, h, Q = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ \eqn{T \times m}{T x m} multiple time
    series object created with \code{ts()}.}
  \item{p}{ Lag length, an integer}
  \item{z}{ NOT IMPLEMENTED AT PRESENT: THIS SHOULD BE A \eqn{T \times
      k}{T x k} matrix of exogenous variables.  Can
    be \code{z = NULL} if there are none (the default). }
  \item{lambda0}{ \eqn{[0,1]}, Overall tightness of the prior
    (discounting of prior scale).}
  \item{lambda1}{ \eqn{[0,1]}, Standard deviation or tightness of the
    prior around the AR(1) parameters. }
  \item{lambda3}{ Lag decay (\eqn{>0}, with 1=harmonic) }
  \item{lambda4}{ Standard deviation or tightness around the intercept \eqn{>0} }
  \item{lambda5}{ Standard deviation or tightness around the exogneous
    variable coefficients \eqn{>0}}
  \item{mu5}{ Sum of coefficients prior weight \eqn{\ge0}.  Larger
    values imply difference stationarity.}
  \item{mu6}{ Dummy initial observations or drift prior \eqn{\ge0}.
    Larger values allow for common trends.}
  \item{nu}{ Prior degrees of freedom, \eqn{m+1} }
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
  \item{prior}{ One of three values: 0 = Normal-Wishart prior, 1 =
    Normal-flat prior, 2 = flat-flat prior (i.e., akin to MLE)}
  \item{h}{ Number of regimes / states, an integer }
  \item{Q}{ \eqn{ h} dimensional transition matrix for the MS
    process. \eqn{h \times h}{h x h} Markov transition matrix whose rows
    sum to 1 with the main weights on the diagonal elements.  Default is
    \code{NULL} and the initial value is defined by \code{qtune}.}
}

\details{ This function sets the initial or starting values for the the
  optimization algorithm for the mode of the MSBVAR models in
  \code{\link{msbvar}}.  This is an attempt to (1) allow for a robust,
  smart guess for starting the block-optimization algorithm and (2)
  allow for user inputs to \code{initialize.opt}.

  The function does three things:

  (1) Estimates an initial \code{\link{szbvar}} model as a baseline,
  non-regime switching model.

  (2) Estimates a set of \code{h} VAR regressions based on a
  \code{\link{kmeans}} clustering of the time series with \eqn{h}
  clusters or centers.  The VAR models fit to each of the \eqn{h}
  subsets of data are used to initialize the \code{\link{msbvar}}
  function. 

  (3) Sets an initial value for \code{Q} in the block optimization
  algorithm for the mode of the MLE / posterior for the MSBVAR model.
  If \code{Q=NULL}, for an \eqn{h \times h}{h x h} transition matrix
  \code{Q}, this initial value is set based on the results from the
  \code{\link{kmeans}} clustering of the data. If the user inputs a
  value of \code{Q}, this is used and error checked to make sure it has
  the correct format (i.e., rows sum to 1, etc.)
  
}
\value{  A list with three elements (these are the inputs for the
  \code{initialize.opt} argument in \code{\link{msbvar}})
  \item{init.model }{ An object of the class SZBVAR, see
    \code{\link{szbvar}} for details}
  \item{thetahat.start }{ The starting values for the regression
    parameters for the block optimization algorithm in
    \code{\link{msbvar}}. This is an \eqn{m \times (mp + 1 + m) \times
      h}{m x (mp+1 + m) x h} array of the initial coefficients.  For the
    \code{i}th element of the array, the
    \eqn{m} rows refer to the equations, the first column elements are
    the intercepts, the next \eqn{2:(mp + 1)} columns are the AR(p)
    coefficients, and the final \eqn{m \times m}{m x m} elements are the
    error covariance for the regime, for that array element.}
  \item{Qhat.start }{ Initial value of \code{Q}}
}
% \references{

% }
\author{
Patrick T. Brandt
}
\note{ 
  This function can be used to model other ways to set the initial
  conditions.  The subsequent calls to the \code{\link{msbvar}} function
  only require an object that satisfies having the elements returned
  from this function --- computed by this function or the user in some
  way. 
}


\seealso{ \code{\link{msbvar}}

}
\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ model }% __ONLY ONE__ keyword per line
