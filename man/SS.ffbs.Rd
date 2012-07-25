\name{SS.ffbs}
\alias{SS.ffbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  State-space forward-filter and backwards-sampler for a
  Markov-switching VAR model}

\description{ This function estimates the \eqn{h} state probabilities for
  all of the observations for a Gaussian likelihood
}
\usage{
SS.ffbs(e, bigt, m, p, h, sig2, Q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{ \eqn{bigt \times m \times h}{TT x m x h} array of the residuals for an
  MSBVAR process}
  \item{bigt}{ integer, number of observations in the model}
  \item{m}{ integer, number of equations or variables in the MSBVAR
    model}
  \item{p}{ integer, number of lags in the model}
  \item{h}{ integer, number of regimes in the MSBVAR model}
  \item{sig2}{ \eqn{m \times m \times h}{m x m x h} array of the
  covariances for each regime (can be the same for each of the \eqn{h} regimes)}
  \item{Q}{ \eqn{h \times h}{h x h} first order Markov transition
  matrix; each row must sum to 1 }
}

\details{
The estimation of an MSBVAR model requires and efficient classifier of
the states for the observed filtered probabilities.  This function
provides a way to accomplish this and is one of the workhorses in the
estimation in the \code{\link{msbvar}} and \code{\link{gibbs.msbvar}}
function.

This function uses compiled Fortran code to draw the 0-1 matrix of the regimes.  It
uses the Baum-Hamilton-Lee-Kim (BHLK) filter and smoother to
estimate the regime probabilities.  Draws are based on the
standard forward-filter-backward-sample algorithm.
}
\value{
  A \eqn{T \times h}{T x h} matrix of the sampled regimes.  Each row
  corresponds to an identity matrix element giving the regime
  classification for the observation.
}
\references{
  Kim, C.J. and C.R. Nelson. 1999.  State-space models with regime
  switching. Cambridge, Mass:  MIT Press.

  Krolzig, Hans-Martin. 1997. Markov-Switching Vector Autoregressions:
  Modeling, Statistical Inference, and Application to Business Cycle
  Analysis.
  
  Sims, Christopher A. and Daniel F. Waggoner and Tao
  Zha. 2008. "Methods for inference in large multiple-equation
  Markov-switching models"  Journal of Econometrics 146(2):255--274.

}
\author{ Patrick T. Brandt }
\note{ This function assumes that the innovation in the MSBVAR model are
  multivariate normal.  The resulting filter and sample follows that in
  the references listed above. This function is provided so users can
  build their own customized
  MSBVAR models.  Users can write functions to generate the (B)VAR
  residuals for their own customized MSBVAR models than then provide the
  residuals and their covariances and transition matrix \eqn{Q}.  This
  function can then be used to estimate / sample the regime
  probabilities. So if you need an MSBVAR model where only certain
  parameters change --- rather than all of them as in the existing
  \code{\link{msbvar}} and \code{\link{gibbs.msbvar}} functions --- you
  can build your own estimator using this function.  This function takes
  care of the hard part of building an MSBVAR model.

}

\seealso{ \code{\link{msbvar}}, \code{\link{gibbs.msbvar}}
}
\examples{
# Simple example to show how data are input to the filter-sampler.
# Assumes a simple bivariate regression model with switching means and
# variances.

TT <- 100
h <- 2
m <- 2
set.seed(123)
x1 <- rnorm(TT)
x2 <- rnorm(TT)
y1 <- 5 + 2*x1 + rnorm(TT)
y2 <- 1 + x2 + 5*rnorm(TT)

Y <- rbind(cbind(y1[1:(0.5*TT)],y2[1:(0.5*TT)]),
           cbind(y2[((0.5*TT)+1):TT],y1[((0.5*TT)+1):TT]))
X <- rbind(cbind(x1[1:(0.5*TT)],x2[1:(0.5*TT)]),
           cbind(x2[((0.5*TT)+1):TT],x1[((0.5*TT)+1):TT]))

u1 <- Y - tcrossprod(cbind(rep(1,TT), X), matrix(c(5,2,0,1,1,0), 2, 3))
u2 <- Y - tcrossprod(cbind(rep(1,TT), X), matrix(c(1,1,0,5,2,0), 2, 3))

u <- array(0, c(TT, m, h))
u[,,1] <- u1
u[,,2] <- u2

Sik <- array(0, c(m,m,h))
Sik[,,1] <- diag(c(1,25))
Sik[,,2] <- diag(c(25,1))

Q <- matrix(c(0.9,0.2,0.1,0.8), h, h)

# estimate the states 100 times
ss <- replicate(100, SS.ffbs(u, TT, m, p=1, h, Sik, Q), simplify=FALSE)

# Get the state estimates from the 100 simulations
ss.est <- matrix(unlist(ss), nrow=(h*TT + h^2))

ss.prob <- matrix(rowMeans(ss.est[1:(h*TT),]), ncol=h)
ss.transition <- matrix(rowMeans(ss.est[((h*TT)+1):((h*TT) + h^2),]),
                        h, h)

}
\keyword{ ts }
\keyword{ models }
