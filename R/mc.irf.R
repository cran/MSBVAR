# Generate the impulse responses -- requires a posterior sample of the
# A0 to be drawn already.  Could probably put an option / overload in here that
# if there are no A0 supplied that it draws them.

"mc.irf" <- function(varobj, nsteps, draws=0, A0.posterior=NULL, sign.list=rep(1,ncol(varobj$Y)))
{
    if(inherits(varobj, "VAR")){
        return(mc.irf.VAR(varobj=varobj, nsteps=nsteps, draws=draws))
    }
    if(inherits(varobj, "BVAR")){
        return(mc.irf.BVAR(varobj=varobj, nsteps=nsteps, draws=draws))
    }
    if(inherits(varobj, "BSVAR")){
        return(mc.irf.BSVAR(varobj=varobj, nsteps=nsteps,
                            A0.posterior=A0.posterior, sign.list=sign.list))
    }
}

"mc.irf.VAR" <- function(varobj, nsteps, draws)
{ .Call("mc.irf.var.cpp", varobj, as.integer(nsteps), as.integer(draws)) }


"mc.irf.BVAR" <- function(varobj, nsteps, draws)
{
    output <- mc.irf.VAR(varobj, nsteps, draws)
    attr(output, "class") <- c("mc.irf.BVAR")
    return(output)
}


"mc.irf.BSVAR" <- function(varobj, nsteps, A0.posterior, sign.list)
{ m<-dim(varobj$ar.coefs)[1]  # Capture the number of variablesbcoefs <- varobj$Bhat
  p<-dim(varobj$ar.coefs)[3]    # Capture the number of lags
  ncoef <- dim(varobj$B.posterior)[1]
  n0 <- varobj$n0
  n0cum <- c(0,cumsum(n0))
  N2 <- A0.posterior$N2

  # Get the covar for the coefficients
  XXinv <- chol(solve(varobj$Hpinv.posterior[[1]]))

  # storage for the impulses and the sampled coefficients.
  impulse <- matrix(0,nrow=N2, ncol=(m^2*nsteps))

  .Call("mc.irf.bsvar.cpp", A0.posterior$A0.posterior,
        as.integer(nsteps), as.integer(N2), as.integer(m),
        as.integer(p), as.integer(ncoef), as.integer(n0),
        as.integer(n0cum), XXinv, varobj$Ui, varobj$P.posterior,
        sign.list)
}

"plot.mc.irf" <- function(x, method=c("Sims-Zha2"), component=1,
                            probs=c(0.16,0.84), varnames=NULL,...)
{
    if(inherits(x, "mc.irf.VAR"))
    { plot.mc.irf.VAR(x, method=method, component=component,
                            probs=probs, varnames=varnames) }

    if(inherits(x, "mc.irf.BVAR"))
    { plot.mc.irf.BVAR(x, method=method, component=component,
                       probs=probs, varnames=varnames) }

    if(inherits(x, "mc.irf.BSVAR"))
    { plot.mc.irf.BSVAR(x, method=method, component=component,
                        probs=probs, varnames=varnames) }

}

"plot.mc.irf.VAR" <- function(x, method, component,
                              probs, varnames=NULL,...)
{ mc.impulse <- x
  m <- sqrt(dim(mc.impulse)[3])
  nsteps <- dim(mc.impulse)[2]
  draws <- dim(mc.impulse)[1]

  # Storage
  irf.ci <- array(0,c(nsteps,length(probs)+1,m^2))

  # Compute the IRF confidence intervals for each element of MxM
  # array

  # There are multiple methods for doing this.

  # Monte Carlo / Bootstrap percentile method
  if (method=="Percentile")
  {
      eigen.sum <- 0
      for(i in 1:m^2)
      {
          irf.bands <- t(apply(mc.impulse[,,i], 2, quantile, probs))
          irf.mean <- apply(mc.impulse[,,i], 2, mean)
          irf.ci[,,i] <- cbind(irf.bands, irf.mean)
      }
  }
  if (method=="Normal Approximation")
  {
      eigen.sum <- 0
      for (i in 1:m^2)
      {
          irf.mean <- apply(mc.impulse[,,i], 2, mean)
          irf.var <- apply(mc.impulse[,,i], 2, var)
          irf.bands <- irf.mean + matrix(rep(qnorm(probs), each=nsteps), nrow=nsteps)*irf.var
          irf.ci[,,i] <- cbind(irf.bands, irf.mean)
      }

  }

  # Sims and Zha symmetric eigen decomposition (assumes normality
  # approximation) with no accounting for the correlation across the
  # responses.  Method does account for the correlation over time.

  if (method=="Sims-Zha1")
  {
      eigen.sum <- matrix(0, m^2, nsteps)
      for(i in 1:m^2)
      {
          decomp <- eigen(var(mc.impulse[,,i]), symmetric=T)
          W <- decomp$vectors
          lambda <- decomp$values
          irf.mean <- apply(mc.impulse[,,i], 2, mean)
          irf.bands <- irf.mean + W[,component]*matrix(rep(qnorm(probs), each=nsteps), nrow=nsteps)*sqrt(lambda[component])
          irf.ci[,,i] <- cbind(irf.bands, irf.mean)
          eigen.sum[i,] <- 100*lambda/sum(lambda)
      }
  }

  # Sims and Zha asymmetric eigen decomposition (no normality
  # assumption) with no accounting for the correlation across the
  # responses.  Method does account for correlation over time.
  if (method=="Sims-Zha2")
  {
      eigen.sum <- matrix(0, m^2, nsteps)
      for(i in 1:m^2)
      {
          decomp <- eigen(var(mc.impulse[,,i]), symmetric=T)
          W <- decomp$vectors
          lambda <- decomp$values
          gammak <- mc.impulse[,,i]*(W[component,])
          gammak.quantiles <- t(apply(gammak, 2, quantile, probs=probs))
          irf.mean <- apply(mc.impulse[,,i], 2, mean)
          irf.bands <- irf.mean + gammak.quantiles
          irf.ci[,,i] <- cbind(irf.bands, irf.mean)
          eigen.sum[i,] <- 100*lambda/sum(lambda)
      }
  }

  # Sims and Zha asymmetric eigen decomposition with no normality
  # assumption and an accounting of the temporal and cross response
  # correlations.

  if (method=="Sims-Zha3")
  {
      eigen.sum <- matrix(0, m^2, m^2*nsteps)

      # Stack all responses and compute one eigen decomposition.
      stacked.irf <- array(mc.impulse, c(draws, m^2*nsteps))
      decomp <- eigen(var(stacked.irf), symmetric=T)
      W <- decomp$vectors
      lambda <- decomp$values
      gammak <- stacked.irf*W[component,]
      gammak.quantiles <- apply(gammak, 2, quantile, probs)
      irf.mean <- matrix(apply(stacked.irf, 2, mean),
                         nrow=length(probs),
                         ncol=dim(stacked.irf)[2],
                         byrow=T)
      irf.bands <- irf.mean + gammak.quantiles

      # Reshape these....
      irf.ci <- array((rbind(irf.bands,irf.mean[1,])),
                      c(length(probs)+1, nsteps, m^2))
      irf.ci <- aperm(irf.ci, c(2, 1, 3))
      eigen.sum <- 100*lambda/sum(lambda)
  }

  # Compute the bounds for the plots
  minmax <- matrix(0, nrow=m, ncol=2)
  within.plots <- apply(irf.ci, 3, range)

  tmp <- (c(1:m^2)%%m)
  tmp[tmp==0] <- m
  indices <- sort(tmp, index.return=T)$ix
  dim(indices) <- c(m, m)

  for(i in 1:m){minmax[i,] <- range(within.plots[,indices[,i]])}
  # Now loop over each m columns to find the minmax for each column
  # responses in the MAR plot.
  j <- 1
  # Plot the results
  par(mfcol=c(m,m),mai=c(0.25,0.25,0.15,0.25), omi=c(0.15,0.75,1,0.15))
  for(i in 1:m^2)
  {
      lims <- ifelse((i-m)%%m==0, m, (i-m)%%m)
      ts.plot(irf.ci[,,i],
              gpars=list(xlab="",ylab="",ylim=minmax[lims,]))
      abline(h=0)

      if(i<=m){ mtext(varnames[i], side=2, line=3)}
      if((i-1)%%m==0){
          mtext(varnames[j], side=3, line=2)
          j <- j+1
      }
  }

  mtext("Response in", side=2, line=3, outer=T)
  mtext("Shock to", side=3, line=3, outer=T)

  # Put response names on the eigenvector fractions
  if(method == "Sims-Zha1" | method == "Sims-Zha2")
  {
      if(is.null(varnames)==T) varnames <- paste("V", seq(1:m), sep = "")
      shock.name <- rep(varnames, m)
      response.name <- rep(varnames, each=m)
      eigen.sum <- cbind(shock.name, response.name, as.data.frame(eigen.sum))
      colnames(eigen.sum) <- c("Shock","Response", paste("Component", seq(1:nsteps)))
  }
  if(method == "Sims-Zha3")
  { names(eigen.sum) <- c(paste("Component", seq(1:m^2*nsteps))) }

  # Invisible return
  invisible(list(responses=irf.ci, eigenvector.fractions=eigen.sum))
}

"plot.mc.irf.BVAR" <- function(x, method, component,
                               probs, varnames, ...)
{
    plot.mc.irf.VAR(x, method, component,
                    probs, varnames=varnames, ...)
}


# BSVAR model IRFs

"plot.mc.irf.BSVAR" <- function(x, method, component, probs, varnames, ...)
{
    m <- sqrt(dim(x)[3])
    nsteps <- dim(x)[2]
    draws <- dim(x)[1]
    irf.ci <- array(0, c(nsteps, length(probs) + 1, m^2))
    if (method == "Percentile") {
        eigen.sum <- 0
        for (i in 1:m^2) {
            irf.bands <- t(apply(x[, , i], 2, quantile,
                probs))
            irf.mean <- apply(x[, , i], 2, mean)
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
        }
    }
    if (method == "Normal Approximation") {
        eigen.sum <- 0
        for (i in 1:m^2) {
            irf.mean <- apply(x[, , i], 2, mean)
            irf.var <- apply(x[, , i], 2, var)
            irf.bands <- irf.mean + matrix(rep(qnorm(probs),
                each = nsteps), nrow = nsteps) * irf.var
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
        }
    }
    if (method == "Sims-Zha1") {
        eigen.sum <- matrix(0, m^2, nsteps)
        for (i in 1:m^2) {
            decomp <- eigen(var(x[, , i]), symmetric = T)
            W <- decomp$vectors
            lambda <- decomp$values
            irf.mean <- apply(x[, , i], 2, mean)
            irf.bands <- irf.mean + W[, component] * matrix(rep(qnorm(probs),
                each = nsteps), nrow = nsteps) * sqrt(lambda[component])
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i, ] <- 100 * lambda/sum(lambda)
        }
    }
    if (method == "Sims-Zha2") {
        eigen.sum <- matrix(0, m^2, nsteps)
        for (i in 1:m^2) {
            decomp <- eigen(var(x[, , i]), symmetric = T)
            W <- decomp$vectors
            lambda <- decomp$values
            gammak <- x[, , i] * (W[component, ])
            gammak.quantiles <- t(apply(gammak, 2, quantile,
                probs = probs))
            irf.mean <- apply(x[, , i], 2, mean)
            irf.bands <- irf.mean + gammak.quantiles
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i, ] <- 100 * lambda/sum(lambda)
        }
    }
    if (method == "Sims-Zha3") {
        eigen.sum <- matrix(0, m^2, m^2 * nsteps)
        stacked.irf <- array(x, c(draws, m^2 * nsteps))
        decomp <- eigen(var(stacked.irf), symmetric = T)
        W <- decomp$vectors
        lambda <- decomp$values
        gammak <- stacked.irf * W[component, ]
        gammak.quantiles <- apply(gammak, 2, quantile, probs)
        irf.mean <- matrix(apply(stacked.irf, 2, mean), nrow = length(probs),
            ncol = dim(stacked.irf)[2], byrow = T)
        irf.bands <- irf.mean + gammak.quantiles
        irf.ci <- array((rbind(irf.bands, irf.mean[1, ])), c(length(probs) +
            1, nsteps, m^2))
        irf.ci <- aperm(irf.ci, c(2, 1, 3))
        eigen.sum <- 100 * lambda/sum(lambda)
    }

    minmax <- matrix(0, nrow = m, ncol = 2)
    within.plots <- apply(irf.ci, 3, range)

    tmp <- (c(1:m^2)%%m)
    tmp[tmp==0] <- m
    indices <- sort(tmp, index.return=T)$ix
    dim(indices) <- c(m, m)

    for(i in 1:m)
      { minmax[i,] <- range(within.plots[,indices[,i]])
      }

    # Now loop over each m columns to find the minmax for each column
    # responses in the MAR plot.
    j <- 1
    # Plot the results
    par(mfcol=c(m,m),mai=c(0.15,0.2,0.15,0.15), omi=c(0.15,0.75,1,0.15))
    for(i in 1:m^2)
      {
          lims <- ifelse((i-m)%%m==0, m, (i-m)%%m)
          ts.plot(irf.ci[,,i],
                  gpars=list(xlab="",ylab="",ylim=minmax[lims,]))

          abline(h=0)

          if(i<=m)
            { mtext(varnames[i], side=2, line=3) }
          if((i-1)%%m==0)
            { mtext(varnames[j], side=3, line=2)
              j <- j+1
            }
        }

    # Add row and column labels for graph

    mtext("Response in", side = 2, line = 3, outer = T)
    mtext("Shock to", side = 3, line = 3, outer = T)


    # Return eigendecomposition pieces if necessary

    if (method == "Sims-Zha1" | method == "Sims-Zha2") {
        if (is.null(varnames) == T)
            varnames <- paste("V", seq(1:m), sep = "")
        shock.name <- rep(varnames, m)
        response.name <- rep(varnames, each = m)
        eigen.sum <- cbind(shock.name, response.name, as.data.frame(eigen.sum))
        colnames(eigen.sum) <- c("Shock", "Response", paste("Component",
            seq(1:nsteps)))
    }
    if (method == "Sims-Zha3") {
        names(eigen.sum) <- c(paste("Component", seq(1:m^2 *
            nsteps)))
    }
    invisible(list(responses = irf.ci, eigenvector.fractions = eigen.sum))
}

