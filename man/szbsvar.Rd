\name{szbsvar}
\alias{szbsvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Structural Sims-Zha Bayesian VAR model estimation}
\description{
  Estimates the posterior mode for a Bayesian Structural Vector
  Autoregression (B-SVAR) model
  using the prior specified by Sims and Zha (1998)
}
\usage{
szbsvar(Y, p, z = NULL,
        lambda0, lambda1, lambda3, lambda4, lambda5,
        mu5, mu6, ident, qm = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ \eqn{T \times m}{T x m} multiple time
    series object created with \code{ts()} with no NAs.}
  \item{p}{ integer lag length for the model }
  \item{z}{ \eqn{T \times k}{T x k} matrix of exogenous variables (not including an
    intercept)}
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
  \item{ident}{ \eqn{m \times m}{m x m} matrix of binary indicators for the identification
    of the free and restricted contemporaneous parameters in \eqn{A_0}{A(0)}. }
  \item{qm}{ Frequency of the data for lag decay equivalence.  Default
    is 4, and a value of 12 will match the lag decay of monthly to
    quarterly data.  Other values have the same effect as "4"}
}
\details{
  This function estimates the posterior mode for the Bayesian structural
  VAR (B-SVAR) model described by Sims and Zha (1998) and Waggoner and
  Zha (2003).  This B-SVAR model is based a specification of the
  dynamic simultaneous equation representation of the model.  The prior
  is constructed for the structural parameters.

  The basic SVAR model has the form of Waggoner and Zha (2003):
  \deqn{ y_t^\prime A_0 = \sum_{\ell=1}^p Y_{t-\ell}^\prime A_\ell +
  z_t^\prime D + \epsilon_t^\prime,  t = 1, \ldots, T,
   }{y(t)' A(0) = \sum_{i=1}^p Y(t-i)' A(i) +
  z(t)' D + e(t)',  t = 1, ..., T,}

  where \eqn{A_i}{A(i)} are \eqn{m \times m}{m x m} parameter matrices for the
  contemporaneous and lagged effects of the endogenous variables, \eqn{D} is
  an \eqn{h \times m}{h x m} parameter matrix for the exogenous variables
  (including an intercept), \eqn{y_t}{y(t)} is the \eqn{m \times 1}{m x
    1} matrix of the endogenous variables, \eqn{z_t}{z(t)} is a \eqn{h
    \times 1}{h x 1} vector of exogenous
  variables (including an intercept) and \eqn{\epsilon_t}{e(t)} is the \eqn{m \times
    1}{m x 1} matrix of structural shocks.  NOTE that in this representation of
  the model, the columns of the \eqn{A_\ell}{A(i)} matrices refer to the equations!

   The structural shocks are normal with mean and variance equal to the
   following:
   \deqn{
   E[\epsilon_t | y_1, \ldots, y_{t-1}, z_1, \ldots z_{t-1}] =
   0}{E[e(t) | Y(1),...,Y(t-1), z(1),...,z(t-1)] = 0}
   \deqn{E[\epsilon_t \epsilon_t^\prime | y_1, \ldots, y_{t-1}, z_1, \ldots
   z_{t-1}] = I}{E[e(t) e(t)' | y(1),..., y(t-1), z(1),...z(t-1)] = I
   }

   The \emph{reduced form representation} of the SVAR model can be found
   by post-multiplying through by \eqn{A_0^{-1}}{A(0)^{-1}}:
   \deqn{
   y_t^\prime A_0 A_0^{-1} = \sum_{\ell=1}^p Y_{t-\ell}^\prime A_\ell
   A_0^{-1} + z_t^\prime DA_0^{-1} + \epsilon_t^\prime
   A_0^{-1}}{y_t' A(0) A(0)^{-1} = \sum_{i=1}^p Y(t-i)' A(i)
   A(0)^{-1} + z_t' D A(0)^{-1} + e(t)' A(0)^{-1}}

   \deqn{y_t^\prime  = \sum_{\ell=1}^p Y_{t-\ell}^\prime B_\ell
   + z_t^\prime \Gamma  + \epsilon_t^\prime A_0^{-1}.}{y(t)'  =
   \sum_{i=1}^p Y(t-i)' B(i) + z(t)' G  + e(t)' A(0)^{-1}.}

  The reduced form error covariance matrix is found from the crossproduct
  of the reduced form innovations:

 \deqn{
  \Sigma = E[(\epsilon_t^\prime A_0^{-1})(\epsilon_t^\prime
  A_0^{-1})^\prime] = [A_0 A_0^\prime]^{-1}.
  }{S = E[(e(t)' A(0)^{-1})(e(t)'A_0^{-1})'] = [A(0) A(0)']^{-1}}.

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

  As in Sims and Zha (1998) and Waggoner and Zha (2003), the prior for
  the model is formed for each of the equations.  To illustrate the
  prior, the model is written in the more compact notation

  \deqn{
  y_t^\prime A_0 = x_t^\prime F + \epsilon_t^\prime}{y(t)' A(0) = x(t)'F
  + e(t)'}
  where
  \deqn{
  x_t^\prime = [ y_{t-1}^\prime \cdots y_{t-p}^\prime, z_t^\prime],
  F^\prime = [A_1^\prime \cdots A_p^\prime \, D^\prime]
  }{x(t)' = [ y(t-1)', ..., y(t-p)', z(t)'],
  F' = [A(1)', ..., A(p)', D']}
  are the matrices of the right hand side variables and the right hand
  side coefficients for the SVAR model.

  The general form of this prior is then
  \deqn{a_i \sim N(0, \bar{S_i}) \quad \textrm{and} \quad f_i | a_i \sim
    N(\bar{P}_i a_i, \bar{H}_i)}{a(i) ~ N(0, S(i)) and f(i) | a(i) ~
    N(P(i) a(i), H(i))}

  where \eqn{\bar{S}_i}{S(i)} is an \eqn{m \times m}{m x m} prior
  covariance of the contemporaneous parameters, and
  \eqn{\bar{H}_i}{H(i)} is the \eqn{k \times k}{k x k}
  prior covariance of the parameters in \eqn{f_i | a_i}{f(i) | a(i)}.
  The prior means of
  \eqn{a_i}{a(i)} are zero in the structural model, while the "random walk"
  component is in \eqn{\bar{P}_i a_i}{P(i) a(i)}.

  The prior covariance matrix of the errors, \eqn{\bar{S}_i}{S(i)}, is
  initially estimated using a VAR(p) model via OLS, with an intercept
  and no demeaning of the data.

  The Bayesian prior is constructed for the unrestricted VAR model and
  then mapped into the restricted prior parameter space, as discussed in
  Waggoner and Zha (2003a).
}

\value{
  A list of the class "BSVAR" that summarizes the posterior mode of the
  B-SVAR model


  \item{XX }{ \eqn{X'X + H_0}{X'X + H_0} crossproduct moment matrix for the predetermined
    variables in the model plus the prior}
  \item{XY}{ \eqn{X'Y} for the model, including the dummy observations for mu5
    and mu6 }
  \item{YY}{ \eqn{m \times m}{m x m} Crossproduct for the Y's in the
    model}
  \item{y}{ \eqn{T \times m}{T x m} input data in \code{dat} plus the
    m dummy observations for \code{dat}}
  \item{structural.innovations}{ \eqn{T \times m}{T x m} structural innovations for the
    SVAR model}
  \item{Ui}{ \eqn{m \times q_i}{m x q(i)} Null space matrices that map
    the columns of \eqn{A_0}{A(0)} to the free parameters of the columns }
  \item{Hpinv.tilde}{ Prior covariance for the predetermined and
    exogenous regression in the B-SVAR}
  \item{H0inv.tilde}{ \eqn{m} dimensional list of the prior covariances for
    the free parameters of the i'th equation in the model's \eqn{A_0} matrix}
  \item{Pi.tilde}{ list of \eqn{(m^2 p + 1 + h) \times q_i}{(m^2 p + 1 +
      h) x q(i)} matrices of the prior for the parameters for the
    predetermined variables in the model }
  \item{Hpinv.posterior}{ \eqn{(m^2 p + 1 + h) \times m}{(m^2 * p + 1 + h) x m} matrix of the
    posterior of the structural parameters for the predetermined variables}
  \item{P.posterior}{ list of \eqn{(m^2 p + 1 + h) \times m}{(m^2 * p +
      1 + h) x m} matrices of the posterior of the paramters for the
    predetermined variables in the model}
  \item{H0inv.posterior}{ \eqn{m} dimensional list of the posterior covariances for
    the free parameters of the i'th equation in the model's \eqn{A_0}{A(0)} matrix}
  \item{A0.mode}{ posterior mode of the \eqn{A_0}{A(0)} matrix}
  \item{F.posterior}{ \eqn{(m^2 p + 1 + h) \times m}{(m^2 * p + 1 + h) x
      m} matrix of the posterior of the structural parameters for the predetermined variables}
  \item{B.posterior}{ \eqn{(m^2 p + 1 + h) \times m}{(m^2 * p + 1 + h) x
      m} matrix of the posterior of the reduced form parameters for the predetermined variables}
  \item{ar.coefs}{ \eqn{(m^2 p) \times m}{m^2 p x m} matrix of the
    posterior of the reduced form autoregressive parameters }
  \item{intercept}{ \eqn{m} dimensional vector of the reduced form
    intercepts }
  \item{exog.coefs}{ \eqn{h \times m}{h x m} matrix of the reduced form exogenous
    variable coefficients }
  \item{prior}{ List of the prior parameter:
    \code{c(lambda0,lambda1,lambda3,lambda4,lambda5, mu5, mu6)}. }
  \item{df}{ Degrees of freedom for the model: T + number of dummy
    observations - lag length. }
  \item{n0}{ \eqn{m} dimensional list of the number of free parameters for the
  \eqn{A_0}{A(0)} matrix for equation \eqn{i}. }
  \item{ident}{ \eqn{m \times m}{m x m} identification matrix
    \code{ident}. }
}
\references{
Sims, C.A. and Tao A. Zha. 1998. "Bayesian Methods for Dynamic
Multivariate Models." \emph{International Economic Review}. 39(4):949-968.

Waggoner, Daniel F. and Tao A. Zha.  2003a. "A Gibbs sampler for
structural vector autoregressions" \emph{Journal of Economic Dynamics
\& Control}. 28:349--366.

Waggoner, Daniel F. and Tao A. Zha. 2003b. "Likelihood preserving
normalization in multiple equation models".  \emph{Journal of
Econometrics}. 114: 329--347.

Brandt, Patrick T. and John R. Freeman. 2006. "Advances in Bayesian
Time Series Modeling and the Study of Politics: Theory Testing,
Forecasting, and Policy Analysis". \emph{Political Analysis} 14(1):1-36.

}
\author{ Patrick T. Brandt }
%\note{ }

\section{Warning }{ If you do not understand the model described here,
  you probably want the models described in \code{\link{szbvar}} or
  \code{\link{reduced.form.var}}}

\seealso{ \code{\link{szbvar}} for reduced form Bayesian VAR models,
  \code{\link{reduced.form.var}} for non-Bayesian reduced form VAR
  models, \code{\link{gibbs.A0}} for drawing from the
  posterior of this model using a Gibbs sampler,
  \code{\link{posterior.fit}} for assessing the
  posterior fit of the model, and \code{\link{mc.irf}}
  for computing impulse responses for this model.
  }
\examples{
# SZ, B-SVAR model for the Levant data
data(BCFdata)
m <- ncol(Y)
ident <- diag(m)
ident[1,] <- 1
ident[2,1] <- 1

# estimate the model's posterior moments
model <- szbsvar(Y, p=2, z=z2, lambda0=0.8, lambda1=0.1,
                 lambda3=1, lambda4=0.1, lambda5=0.05,
                 mu5=0, mu6=5, ident, qm=12)

}
\keyword{ ts }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
