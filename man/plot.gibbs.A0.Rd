\name{plot.gibbs.A0}
\alias{plot.gibbs.A0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a parameter density summary for B-SVAR A(0) objects }
\description{
  Generates an \eqn{m \times m}{m x m} matrix of density plots for each
  free parameter in an \code{szbsvar} \eqn{A_0}{A(0)} object produced by
  \code{gibbs.A0}, with associated highest posterior density (HPD) regions.
}
\usage{
\method{plot}{gibbs.A0}(x, hpd = 0.68, varnames=attr(x, "eqnames"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An \eqn{A_0}{A(0)} posterior object created by
    \code{szbsvar}. }
  \item{hpd}{ Probability width of the highest posterior density region,
  default is 0.68 or approximately one standard deviation around the
  mode of the parameter}
  \item{varnames}{ List of variable names for labeling the equations
      and variables.  Default are the names of the variables for the
      input data to \code{\link{szbsvar}} as fed through
      \code{\link{gibbs.A0}}.  For an SVAR, users often want to relabel
      these as economic sectors or groups of actors for the time series
      and this is the place this can be done.} 
  \item{\dots}{ optional graphics arguments }
}

\details{
  This function plots an \eqn{m \times m}{m x m} matrix of densities for
  the posterior of the \eqn{A_0}{A(0)} free parameters for a B-SVAR
  model.  The plot is arranged such that the unrestricted parameters for
  each contemporaneous effect of each variable on an equation are in the
  row for that equation.  So the first row shows densities for the
  contemporaneous effects of the column variables (as in an impulse
  response plot like \code{plot.irf} or \code{plot.mc.irf}).  Elements
  of  \eqn{A_0}{A(0)} that were restricted to zero are left empty in the
  matrix of densities.  The pattern of the densities will match the
  *tranpose* of the \code{ident} matrix passed to \code{szbsvar}.

  Highest posterior density regions are plotted using Hyndman's 91996)
  density quantile algorithm.  These HPDs are defined by a set of
  vertical bars over the HPD interval.  The vertical line in each plot
  measures the value of the density at the boundaries of the HPD
  region. The HDR is superimposed at the bottom of each density.
}
\value{
  None.  Main purposed is to plot density summaries and HPDs for each of
  the free parameters in an \eqn{A_0}{A(0)} matrix.
}
\references{ Hyndman, Rob J. 1996. "Computing and Graphic Highest
  Density Regions", The American Statistician, 50(2):120--126

  HPD code is borrowed from Hyndman's hdrcde package, version 2.07.
 }
\author{ Patrick T. Brandt }
\note{ The plot will tend to be large, so be sure to adjust the size of
 your plotting device accordingly so things are visible.}
\seealso{ \code{\link[coda]{plot.mcmc}},
  \code{\link[coda]{summary.mcmc}}, and \code{\link{A02mcmc}}.
}  
\examples{
# SZ, B-SVAR model for the Levant data
data(BCFdata)
m <- ncol(Y)
ident <- diag(m)
ident[1,] <- 1
ident[2,1] <- 1

# estimate the model's posterior moments
set.seed(123)
model <- szbsvar(Y, p=2, z=z2, lambda0=0.8, lambda1=0.1, lambda3=1, lambda4=0.1,
                 lambda5=0.05, mu5=0, mu6=5, ident, qm=12)

# Set length of burn-in and size of posterior.  These are only an
# example.  Production runs should set these much higher.
N1 <- 1000
N2 <- 10000

A0.posterior.obj <- gibbs.A0(model, N1, N2, thin=1)

# Plot the matrix of the densities
dev.new()
plot.gibbs.A0(A0.posterior.obj, hpd=0.68, varnames=colnames(Y))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ manip }% __ONLY ONE__ keyword per line
