\name{simulateMSVAR}
\alias{simulateMSVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a Markov-switching VAR (MSVAR) process
}

\description{  Simulate Markov-switching vector autoregression data
}
\usage{
simulateMSVAR(bigt, m, p, var.beta0, var.betas, e.vcv, Q, seed = 214)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bigt}{ Integer, number of observations to generate.
}
  \item{m}{ Integer, number of equations in the VAR process
}
  \item{p}{ Integer, lag length of the VAR(p) process.
}
  \item{var.beta0}{ Array of dimension \eqn{m \times 1 \times h}{m x 1 x
      h} of the VAR intercepts for each regime (h)
}
  \item{var.betas}{ Array of dimension \eqn{m \times mp \times h}{m x mp
      x h} of the autoregressive coefficients.  In each element of the
    array, rows correspond to
    equations, columns to lags.  The first \eqn{m \times m}{m x m}
    columns are the AR(1) coefficients, etc.  
}
  \item{e.vcv}{ Array of dimentsion \eqn{m \times m \times h}{ m x m x
      h} of the error covariances.  The \eqn{m \times m}{  m x m}
    matrices are the error covariances for each regime.
}
  \item{Q}{ \eqn{ h} dimensional transition matrix for the MS
    process. \eqn{h \times h}{h x h} Markov transition matrix whose rows
    sum to 1 with the main weights on the diagonal elements.
}
  \item{seed}{  Integer.  Random number seed.
}
}
\details{ This function simulates a multivariate Markov-switching model,
  MSVAR with \eqn{m} equations, \eqn{p} lags and \eqn{h} regimes. The
  assumption is that the error process is Gaussian.
}
\value{
  A list with two elements:
  \item{Y }{ The simulated MSVAR time series in a \code{ts} object of
    dimension \eqn{ bigt \times m}{ bigt x m}.}
  \item{st }{ A vector of integers identifying the regime of each
    observation in \code{Y}}
}
\references{
  Kim, Chang-Jin and Charles R. Nelson. 1999. State-Space Models with
  Regime Switching: Classical and Gibbs-Sampling Approaches with
  Applications. Cambridge: MIT Press.

}
\author{
  Patrick T. Brandt and Ryan Davis
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{simulateMSAR}} for the univariate version;
  \code{\link{msbvar}}
}

\examples{
\dontrun{
# Example: simulate an MS(h)-VAR(p) model with two equations.
# Have h = 2, m=2, and p=1, simplest case

# VAR simulation parameters
bigt <- 500  # number of observations
m <- 2  # number of endogenous variables
p <- 1  # lag length
h <- 2  # number of regimes

# setup transition matrix with two states

Q <- matrix(c(.98, .02,
              .05, .95), nrow=h, byrow=TRUE)

# theta stores paramater values
# 1st column is intercept
# 2:m*p are the AR coefficients
# (mp+2)'th columns are variance

# regime 1
var.beta0.st1  <- c(1,2)  # intercepts
var.betas.st1  <- matrix(c(.7, .1,
                           .1, .7), m, byrow=TRUE)
# regime 2
var.beta0.st2  <- c(0,0)  # intercepts
var.betas.st2  <- matrix(c(.2, .1,
                           .2, .1), m, byrow=TRUE)

# Build the array
var.beta0 <- array(NA, c(m,1,h))
var.betas <- array(NA, c(m,p*m,h))
var.beta0[,,1] <- var.beta0.st1
var.beta0[,,2] <- var.beta0.st2
var.betas[,,1] <- var.betas.st1
var.betas[,,2] <- var.betas.st2

# Variance-Covariance Matrix for
# MVN distributed disturbances
# regime 1
e.vcv.st1  <- matrix(c(.3, .1,
                       .1, .3), 2)
# regime 2
e.vcv.st2  <- matrix(c(.1, .05,
                       .05, .1), 2)
# combine
e.vcv <- array(NA, c(m, m, h))
e.vcv[,,1] <- e.vcv.st1
e.vcv[,,2] <- e.vcv.st2

# hold true values of parameters for easy comparison to estimates
theta.true.var <- array(NA, c(m, 1+m*p+m, h))
theta.true.var[,1,] <- var.beta0
theta.true.var[,2:(1+p*m),] <- var.betas
theta.true.var[,(1+m*p+1):ncol(theta.true.var),] <- e.vcv

simdata <- simulateMSVAR(bigt, m, p, var.beta0, var.betas, e.vcv, Q)

# Plot
plot(as.ts(simdata[[1]]))

# Fit a simple model

model <- msvar(Y=simdata[[1]], p=1, h=2, niterblkopt=50)

# Plot regime estimates and compare to true simulated values
par(mfrow=c(2,1))
plot(ts(model$fp))
plot(ts(simdata$st))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }
\keyword{ model }% __ONLY ONE__ keyword per line
