\name{gibbs.msbvar}
\alias{gibbs.msbvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gibbs sampler for a Markov-switching Bayesian reduced form
  vector autoregression model }
\description{
  Draws a Bayesian posterior sample for a Markov-switching Bayesian
  reduced form vector autoregression model based on the setup from the
  \code{\link{msbvar}} function.
}
\usage{
gibbs.msbvar(x, N1 = 1000, N2 = 1000, permute = TRUE,
             Beta.idx = NULL, Sigma.idx = NULL, Q.method="MH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ MSBVAR setup and posterior mode estimate generated using the
    \code{\link{msbvar}} function.}
  \item{N1}{ Number of burn-in iterations for the Gibbs sampler (should
    probably be greater than or equal to 1000)}
  \item{N2}{ Number of iterations in the posterior sample. }
  \item{permute}{ Logical (default = TRUE).  Should random permutation
    sampling be used to explore the h! posterior modes? }
  \item{Beta.idx}{ A two element vector indicating the MSBVAR
    ceofficient matrix that is to be ordered for non-permutation
    sampling, i.e., the ordering of the states.  The states will be put
    into ascending order for the parameter selected.  The two elements
    provide are for the two-dimensional array of the VAR
    coefficients. The first number gives the coefficient, the second the
    equation numbers.  Coefficients are ordered by lag, then variable.
    So for an \eqn{m} equation VAR where we want the AR(1) coefficient on
    the second variable's equation, use \code{c(2,2)}.  The intercept is
    the last value, or \eqn{mp+1}.  So the intercept for the first
    equation in a 4 variable model with two lags is \code{c(9,1)}.}
  \item{Sigma.idx}{Scalar integer giving the equation variance that is
    to be ordered for non-permutation sampling, i.e., the ordering of
    the states. The states will be put
    into ascending order for the variance parameter selected.  So if you
    want to identify the results based on equation three, set
    \code{Sigma.idx=3}}
  \item{Q.method}{ choice of the sampler step for the transition matrix,
    \code{Q}. \code{default=MH} uses a Metropolis-Hastings algorithm
    that assumes a stationary Markov process.  The other option is
    \code{Gibbs} which uses a Gibbs sampler Dirichlet draw for a
    non-stationary Markov-switching process.  See Fruwirth-Schnatter
    (2006: 318, 340-341 for details)}
}

\details{
  This function implements a Gibbs sampler for the posterior of a MSBVAR
  model setup with \code{\link{msbvar}}.  This is a reduced form MSBVAR
  model.  The estimation is done in a mixture of native R code and
  Fortran. The sampling of the BVAR coefficients, the transition matrix,
  and the error covariances for each regime are done in native R code.
  The forward-filtering-backward-sampling of the Markov-switching process
  (The most computationally intensive part of the estimation) is handled
  in compiled Fortran code.  As such, this model is reasonably fast for
  small samples / small numbers of regimes (say less than 5000
  observations and 2-4 regimes).  The reason for this mixed
  implementation is that it is easier to setup variants of the model
  (E.g., Some coefficients switching, others not; different sampling
  methods, etc.  Details will come in future versions of the package.)

  The random permuation of the states is done using a multinomial step:
  at each draw of the Gibbs sampler, the states are permuted using a
  multinomial draw.  This generates a posterior sample where the states
  are unidentified.  This makes sense, since the user may have little
  idea of how to select among the \eqn{h!} posterior models of the
  reduced
  form MSBVAR model (see e.g., Fruhwirth-Schnatter (2006)).  Once a
  posterior sample has been draw with random permuation, a clustering
  algorithm (see \code{\link{plotregimeid}}) can be used to identify the
  states, for example, by
  examining the intercepts or covariances across the regimes (see the
  example below for details).

  Only the \code{Beta.idx} or \code{Sigma.idx} value is followed.  If
  the first is given the second will be ignored.  So variance ordering
  for identification can only be used when \code{Beta.idx=NULL}.  See
  \code{\link{plotregimeid}} for plotting and summary methods for the
  permuted sampler.

  The Gibbs sampler is estimated using six steps:
  \describe{
    \item{Drawing the state-space for the Markov process }{ This step
      uses compiled code to draw the 0-1 matrix of the regimes.  It
      uses the Baum-Hamilton-Lee-Kim (BHLK) filter and smoother to
      estimate the regime probabilities.  Draws are based on the
      standard forward-filter-backward-sample algorithm.}
    \item{Drawing the Markov transition matrix \eqn{Q} }{ Conditional on
      the other parameters, this takes a draw from a Dirichlet posterior
      with the \code{alpha.prior} prior.}
    \item{Regression step update }{Conditional on the state-space and
      the Markov-switching process data augmentation steps, estimate a
      set of \eqn{h} regressions, one for each regime.}
    \item{Draw the error covariances, \eqn{\Sigma_h}{Sigma(h)} }{
      Conditional on the other steps, compute and draw the error
      covariances from an inverse Wishart pdf.}
    \item{Draw the regression coefficients }{For each set of classified
      observations' (based on the previous step) BVAR regression
      coefficients, take a draw from their multivariate normal
      posterior.}
    \item{Permute the states }{If \code{permute = TRUE}, then permute
      the states and the respective coefficients.}
  }
  
  The state-space for the MS process is a \eqn{T \times h}{T x h} matrix
  of zeros and ones.  Since this matrix classifies the observations
  infor states for the \code{N2} posterior draws, it does not make sense
  to store it in double precisions.  We use the \code{\link[bit]{bit}}
  package to compress this matrix into a 2-bit integer representation
  for more efficient storage.  Functions are provided (see below) for
  summarizing and plotting the resulting state-space of the MS process.

%  Talk about permutation and why we need it and why you need to estimate
%  the model twice!
}  
\value{
  A list summarizing the reduced form MSBVAR posterior:
  \item{Beta.sample }{ \eqn{N2 \times h(m^2 p + m)}{N2 x h(m^2 p + m)}
    of the BVAR regression coefficients for each regime.  The ordering
    is based on regime, equation, intercept (and in the future
    covariates).  So the first \eqn{p} coefficients are the the first
    equation in the first regime, ordered by lag, not variable; the next
    is the intercept.  This pattern repeats for the remaining
    coefficents across the regimes.}
  \item{Sigma.sample }{\eqn{N2 \times h(\frac{m(m+1)}{2}) }{N2 x
      0.5*h(m(m+1))} matrix of the covariance parameters for the error
    covariances \eqn{\Sigma_h}{\Sigma(h)}.  Since these matrices are
    symmetric p.d., we only store the upper (or lower) portion.  The
    elements in the matrix are the first, second, etc. columns / rows of
    the lower / upper version of the matrix.}
  \item{Q.sample }{ \eqn{N2 \times h^2}{N2 x h^2} }
  \item{transition.sample }{ An array of \code{N2} \eqn{h \times h}{h x
    h} transition matrices.}
  \item{ss.sample }{ List of class SS for the \code{N2} estimates of the state-space
    matrices coded as \code{\link[bit]{bit}} objects for compression /
    efficiency.}
  \item{pfit}{ A list of the posterior fit statistics for the MSBVAR
    model. }
  \item{init.model}{ Initial model -- a varobj from a BVAR like
    \code{\link{szbvar}} that sets up the data and priors.  See
    \code{\link{szbvar}} for a description.}
  \item{alpha.prior}{ Prior for the state-space transitions Q.  This is
    set in the call to \code{\link{msbvar}} and inherited here.}
  \item{h}{ integer, number of regimes fit in the model.}
  \item{p}{ integer, lag length}
  \item{m}{ integer, number of equations}
}
\references{
  Brandt, Patrick T. 2009. "Empirical, Regime-Specific Models
  of International, Inter-group Conflict, and Politics"
 
  Fruhwirth-Schnatter, Sylvia. 2001. "Markov Chain Monte Carlo
  Estimation of Classical and Dynamic
  Switching and Mixture Models". Journal of the American Statistical
  Association. 96(153):194--209.

  Fruhwirth-Schnatter, Sylvia. 2006. Finite Mixture and Markov Switching
  Models. Springer Series in Statistics New York: Springer.
 
  Sims, Christopher A. and Daniel F. Waggoner and Tao
  Zha. 2008. "Methods for inference in large multiple-equation
  Markov-switching models"  Journal of Econometrics 146(2):255--274.
  
  Krolzig, Hans-Martin. 1997. Markov-Switching Vector Autoregressions:
  Modeling, Statistical Inference, and Application to Business Cycle
  Analysis. 
  }
\author{ Patrick T. Brandt }
\note{

  Users need to call this function twice (unless they have really good
  a priori identification information!)  The first call will be using
  the random permutation sampler (so with \code{permute = TRUE}) and
  then some exploration of the clustering of the posterior.  Then, once
  the posterior is identified (i.e., you have chosen one of the \eqn{h!}
  posterior modes), the function is called with \code{permute = FALSE}
  and values specified for \code{Beta.idx} or \code{Sigma.idx}.  See the
  example below for usage.
  
}
\seealso{ \code{\link{msbvar}} for initial mode finding,
  \code{\link{plot.SS}} for plotting regime probabilities,
  \code{\link{mean.SS}} for computing the mean regime probabilities,
  \code{\link{plotregimeid}} for identifying the regimes from a permuted
  sample. }

\examples{
\dontrun{
# This example can be pasted into a script or copied into R to run.  It
# takes a few minutes, but illustrates how the code can be used

data(IsraelPalestineConflict)  

# Find the mode of an msbvar model
# Initial guess is based on random draw, so set seed.
set.seed(123)

xm <- msbvar(IsraelPalestineConflict, p=3, h=2,
             lambda0=0.8, lambda1=0.15,
             lambda3=1, lambda4=1, lambda5=0, mu5=0,
             mu6=0, qm=12,
             alpha.prior=matrix(c(10,5,5,9), 2, 2))

# Plot out the initial mode
plot(ts(xm$fp))
print(xm$Q)

# Now sample the posterior
N1 <- 1000
N2 <- 2000

# First, so this with random permutation sampling
x1 <- gibbs.msbvar(xm, N1=N1, N2=N2, permute=TRUE)

# Identify the regimes using clustering in plotregimeid()
plotregimeid(x1, type="all")

# Now re-estimate based on desired regime identification seen in the
# plots.  Here we are using the intercept of the first equation, so
# Beta.idx=c(7,1).

x2 <- gibbs.msbvar(xm, N1=N1, N2=N2, permute=FALSE, Beta.idx=c(7,1))

# Plot regimes
plot.SS(x2)

# Summary of transition matrix
summary(x2$Q.sample)

# Plot of the variance elements
plot(x2$Sigma.sample)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ts }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
