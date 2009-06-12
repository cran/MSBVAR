\name{msbsvar}
\alias{msbsvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Markov-Switching Sims-Zha Bayesian VAR Model estimation}
\description{
  Estimates the posterior model for a Markov-Switching
  Bayesian Structural Vector Autoregression (B-SVAR) model
  using the prior specified by Sims and Zha (1998)
}
\usage{
msbsvar(Y, z = NULL, p, h, ident,
        lambda0, lambda1, lambda3, lambda4, lambda5, mu5, mu6,
        qm, alpha.prior, max.iter = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ \eqn{T \times m}{T x m} multiple time
    series object created with \code{ts()} with no NAs.}
  \item{z}{ \eqn{T \times k}{T x k} matrix of exogenous variables (not including an
    intercept) CURRENTLY not implemented}
  \item{p}{ integer lag length for the model }
  \item{h}{ \eqn{h} = number of states of regimes in the
    Markov-switching (MS) process for the model parameters }
  \item{ident}{ \eqn{m \times m}{m x m} matrix of binary indicators for the identification
    of the free and restricted contemporaneous parameters in \eqn{A_0(s_t)}{A(0,s(t))}. }
  \item{lambda0}{ \eqn{[0,1]}, Overall tightness of the prior (discounting
    of prior scale). }
  \item{lambda1}{ \eqn{[0,1]}, Standard deviation or tightness of the
    prior around the AR(1) parameters. }
  \item{lambda3}{ Lag decay (\eqn{> 0}, with 1=harmonic) }
  \item{lambda4}{ Standard deviation or tightness around the intercept \eqn{>0} }
  \item{lambda5}{ Standard deviation or tightness around the exogneous
    variable coefficients \eqn{>0}}
  \item{mu5}{ Sum of coefficients prior weight \eqn{\ge0}.  Larger
    values imply difference stationarity.}
  \item{mu6}{ Dummy Initial observations or drift prior \eqn{\ge
      0}.  Larger values allow for common trends.}
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
  \item{alpha.prior}{ \eqn{h \times h}{h x h} matrix for the prior for
    the Dirichlet process for the Markov process matrix \eqn{Q} (see below).}
  \item{max.iter}{ Maximum number of outer iterations for the block
    optimizer for the model parameters.  Default is 10. }
}
\details{
  This function estimates the posterior mode for a version of the Markov-switching
  Bayesian structural VAR (MSBSVAR) model described by Sims, Waggoner,
  and Zha (2008).  This MSBSVAR model is based on a specification of the
  dynamic simultaneous equation representation of the model.  The prior
  is constructed for the structural parameters and an unrestricted
  Markov process (contra the main results in Sims et al. (2008).

  The basic MSBSVAR model has the form of Waggoner and Zha (2003) and
  Sims et al (2008):
  \deqn{ y_t^\prime A_0(s_t) = \sum_{\ell=1}^p Y_{t-\ell}^\prime A_\ell(s_t) +
  z_t^\prime D(s_t) + \epsilon_t(s_t)^\prime,  t = 1, \ldots, T,
   }{y(t)' A(0, s(t)) = \sum_{i=1}^p Y(t-i)' A(i,s(t)) +
  z(t)' D(s(t)) + e(t, s(t))',  t = 1, ..., T,}

  where \eqn{A_i(s_t)}{A(i, s(t))} are \eqn{m \times m}{m x m} parameter matrices for the
  contemporaneous and lagged effects of the endogenous variables in
  regime \eqn{s_t}{s(t)}, \eqn{D(s_t)}{s(t)} is
  an \eqn{h \times m}{h x m} parameter matrix for the exogenous variables
  (including an intercept) in regime \eqn{s_t}{s(t)}, \eqn{y_t}{y(t)} is the \eqn{m \times 1}{m x
    1} matrix of the endogenous variables, \eqn{z_t}{z(t)} is a \eqn{h
    \times 1}{h x 1} vector of exogenous
  variables (including an intercept) and \eqn{\epsilon_t(s_t)}{e(t, s(t))} is the \eqn{m \times
    1}{m x 1} matrix of structural shocks in regime \eqn{s_t}{s(t)}.  NOTE that in this representation of
  the model, the columns of the \eqn{A_\ell(s_t)}{A(i, s(t))} matrices refer to the equations!
   
   The structural shocks are normal with mean and variance equal to the
   following:
   \deqn{
   E[\epsilon_t | y_1, \ldots, y_{t-1}, z_1, \ldots z_{t-1}] =
   0}{E[e(t) | Y(1),...,Y(t-1), z(1),...,z(t-1)] = 0}
   \deqn{E[\epsilon_t \epsilon_t^\prime | y_1, \ldots, y_{t-1}, z_1, \ldots
   z_{t-1}] = I}{E[e(t) e(t)' | y(1),..., y(t-1), z(1),...z(t-1)] = I
   }
   
   At present, the model does NOT include switching structural variances
   (contra the Sims et al. (2008) paper).  This is because this depends
   on the normalization of the variances across the regimes and we are
   working on a more general specification of this for users.

   Restrictions on the contemporaneous parameters in \eqn{A_0}{A(0)} are
   expressed by the specification of the \code{ident} matrix that defines
   the shocks that "hit" each equation in the contemporaneous
   specification.  If \code{ident} is defined as in the following table,

   \tabular{lccc}{
     \tab       \tab Equations \tab \cr
     Variables \tab Eqn 1 \tab Eqn 2 \tab Eqn 3\cr
     Var. 1    \tab 1 \tab  0 \tab  0\cr
     Var. 2    \tab 1 \tab  1 \tab  0\cr
     Var. 3    \tab 0 \tab  1 \tab  1
   }

   then the corresponding \eqn{A_0}{A(0)} is restricted to
   \tabular{lccc}{
     \tab       \tab Equations \tab \cr
     Variables \tab Eqn 1 \tab Eqn 2 \tab Eqn 3\cr
     Var. 1 \tab \eqn{a_{11}}{a(11)} \tab  0 \tab  0\cr
     Var. 2 \tab \eqn{a_{12}}{a(12)} \tab  \eqn{a_{22}}{a(22)} \tab  0\cr
     Var. 3 \tab 0 \tab  \eqn{a_{23}}{a(23)} \tab  \eqn{a_{33}}{a(33)}
   }

  which is interpreted as shocks in variables 1 and 2 hit equation 1
  (the first column); shocks in variables 2 and 3 hit the second
  equation (column 2); and, shocks in variable 3 hit the third equation
  (column 3).

  Note that the identification is the same across the regimes.  (I hope
  to allow different structures across regimes in future versions).
  
  The prior needs to be described here -- borrow from the earlier MPSA
  paper!
  
}
\value{
  A list of the class "MSBSVAR" that summarizes the posterior mode of the MSBSVAR model
  
  \item{b}{ Reduced rank vector of parameters that describe the elements
    of \eqn{A_0}{A(0)}.}
  \item{xi}{ \eqn{m \times h}{m x h} matrix of the regime dependent
    variance weights. }
  \item{Q}{ \eqn{h \times h}{h x h} matrix estimate of the Markov transition
    matrix}
  \item{fp}{ \eqn{T \times h}{T x h} matrix of the filter probabilities
    for the h regimes.}
  \item{m}{ \eqn{m} number of endogenous variables in the system}
  \item{p}{ \eqn{p} number of lags used in the VAR}
  \item{h}{ \eqn{h} number of regimes}
  \item{init.model}{ an object of the class "BSVAR" that contains the
    information for the model setup, prior, etc. (needed for later calculations)}
  \item{n0}{ number of free parameters in each column of the
    \eqn{A_0(s_t)}{A(0, s(t))} matrices}.
  \item{n0cum}{ a cumulative sum of n0}
  \item{abar}{ \eqn{m \times h}{m x h} matrix of the alpha parameters
    for the gamma prior for regime dependent variance weights. }
  \item{bbar}{  \eqn{m \times h}{m x h} matrix of the beta parameters
    for the gamma prior for regime dependent variance weights. }
  \item{alpha.prior}{\eqn{h \times h}{h x h} matrix of values for the
    Dirichlet prior for the Markov transition matrix}

}
\references{
  Sims, C.A., D. Waggoner and T. A. Zha. 2008. "Methods for Inference in
  Large Multiple-Equation Markov-switching Models." \emph{Journal of
  Econometrics}. 
  
  Sims, C.A. and Tao A. Zha. 1998. "Bayesian Methods for Dynamic
  Multivariate Models." \emph{International Economic Review}. 39(4):949-968.

  Waggoner, Daniel F. and Tao A. Zha.  2003a. "A Gibbs sampler for
  structural vector autoregressions" \emph{Journal of Economic Dynamics
  \& Control}. 28:349--366.

Waggoner, Daniel F. and Tao A. Zha. 2003b. "Likelihood preserving
normalization in multiple equation models".  \emph{Journal of
Econometrics}. 114: 329--347.

}
\author{ Patrick T. Brandt}
%\note{ }
\seealso{ \code{\link{gibbs.msbsvar}} for sampling from the model's
  posterior, \code{\link{szbsvar}} for a non-Markov-switching version of
  this model, \code{\link{szbvar}} for reduced form Bayesian VAR models,
  \code{\link{reduced.form.var}} for non-Bayesian reduced form VAR
  models }

% \examples{}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ ts }
\keyword{ regression }% __ONLY ONE__ keyword per line
