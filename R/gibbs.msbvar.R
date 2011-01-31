# gibbs.msbvar.R Functions for Gibbs sampling an MSBVAR model based on
# a SZ prior.

# Patrick T. Brandt
# 20081113 : Initial version
# 20100325 : Added the computation for the log marginal data densities
# 20100615 : Cleaned up permutation / flipping code.  This now has its
#            own function that is called inside the Gibbs sampler.
#            Also reorganized the code in this file so that is easier
#            to read, document and follow.
# 20100617 : Wrote separate functions for the Beta, Sigma, and e block
#            updates.  Cleans up the code for the gibbs.msbvar into
#            managable chunks.



####################################################################
# Utility functions for managing matrices in the Gibbs sampler
####################################################################
#
# vech function for efficiently sorting and subsetting the unique
# elements of symmetric matrices (like covariance Sigma)
#
vech <- function (x)
{
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
        stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
}

####################################################################
# xpnd function for reconstructing symmetric matrices from vech
####################################################################

xpnd <- function (x, nrow = NULL)
{
    dim(x) <- NULL
    if (is.null(nrow))
        nrow <- (-1 + sqrt(1 + 8 * length(x)))/2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag = TRUE)] <- 0
    output <- output + t(hold)
    return(output)
}

####################################################################
# Sampler blocks for the MSBVAR Gibbs sampler
####################################################################
# Metropolis-Hastings sampler for the transition matrix Q
MSBVAR.Q.draw <- function(Q, SStrans, prior)
{
    Q.old <- Q
    h <- dim(Q)[1]
    # compute the density ordinate based on the counts.
    alpha <- SStrans + prior - matrix(1, h, h)

    Q.new <- t(apply(alpha, 1, rdirichlet, n=1))

    # now form the components for the MH ratio
    pnew <- ddirichlet(Q.new, SStrans + matrix(1, h, h))
    pold <- ddirichlet(Q.old, SStrans + matrix(1, h, h))
    gnew <- ddirichlet(Q.new, alpha)
    gold <- ddirichlet(Q.old, alpha)

    A <- sum(log(pnew, gold)) - sum(log(pold, gnew))
    A <- ifelse(is.na(A)==TRUE, -100, A)

    if( log(runif(1)) > A ) return(Q.old) else return(Q.new)
}

####################################################################
# State-space sampling for the regimes -- this is the FFBS
# implementation.  Actual work is done by calling compiled C code in
# the package
####################################################################

SS.draw <- function(u, TT, m, h, Sik, Q)
{
    u.vec <- matrix(0, (TT*m), h); sig.vec <- matrix(0,m*m,h)
    for(i in 1:h){
        u.vec[,i] <- as.vector(t(u[,,i]))
        sig.vec[,i] <- as.vector(Sik[,,i])
    }

    # Compute the steady state / ergodic Q for the initialization.
    SQ <- steady.Q(Q)

    # Sample the state space sampler
    ss <- .Call("SS.newdraw.cpp", u.vec, sig.vec, Q, SQ,
                as.integer(TT), as.integer(h), as.integer(m))
    return(ss)
}

####################################################################
# Sample the error covariances for each state
####################################################################

Sigma.draw <- function(m, h, ss, hreg, Sigmai)
{
    # Draw the variances for the MSVAR from an inverse Wishart --
    # draws h matrices which are returned as an m^2 x h matrix.  These
    # then need to be unwound into the matrices we need.

    # Unscaled inverse Wishart draws
    df <- colSums(ss$SS)
    tmp <- array(unlist(lapply(sapply(1:h,
                                      function(i) {rwishart(N=1,
                                                            df=df[i]+m,
                                                            Sigma=diag(m))},
                                      simplify=FALSE), solve)), c(m,m,h))

    # Scale the draws for each regime
    for (i in 1:h)
    {
        Sigmai[,,i] <- t(chol(hreg$Sigma[,,i]))%*%(tmp[,,i]*(df[i]))%*%chol(hreg$Sigma[,,i])
    }

    return(Sigmai)
}

####################################################################
# Sample the regressors for each state
####################################################################

Beta.draw <- function(m, p, h, Sigmai, hreg, init.model, Betai)
{
    for(i in 1:h) {
        # Get the variance of the regressors
        bcoefs.se <- t(chol(kronecker(Sigmai[,,i],
                                      solve(hreg$moment[[i]]$Sxx +
                                            init.model$hstar))))
        # Sample the regressors
        Betai[,,i] <- hreg$Bk[,,i] + matrix(bcoefs.se%*%matrix(rnorm(m^2*p+m), ncol=1), ncol=m)
    }

    return(Betai)
}

####################################################################
# Update the residuals for each state
####################################################################
residual.update <- function(m, h, init.model, Betai, e)
{
    for(i in 1:h)
    {
        e[,,i] <- init.model$Y[(m+1):nrow(init.model$Y),] - init.model$X[(m+1):nrow(init.model$X),]%*%Betai[,,i]
    }
    return(e)
}

######################################################################
# Function to do the random permuation / state labeling steps for the
# models.  This function handles both the random permutation steps and
# the identified models based on either tbe Beta.idx or the Sigma.idx
# input args.
######################################################################
PermuteFlip <- function(x, h, permute, Beta.idx, Sigma.idx)
{
    # Random permutation step
    if(permute==TRUE){
        rp <- sort(runif(h), index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai; Betaitmp <- x$Betai
        etmp <- x$e; sstmp <- x$ss

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,rp[i]]
            Betaitmp[,,i] <- x$Betai[,,rp[i]]
            etmp[,,i] <- x$e[,,rp[i]]
            sstmp$SS[,i] <- x$ss$SS[,rp[i]]
        }

        Q <- x$Q[rp,rp]
        sstmp$transitions <- x$ss$transitions[rp, rp]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp, ss=sstmp, e=etmp, Q=Q))
    }

    # Sort regimes based on values of regressors in Beta using the
    # Beta.idx indices

    if(permute==FALSE & is.null(Beta.idx)==FALSE){
        # Find the sort indices for the objects
        ix <- sort(x$Betai[Beta.idx[1], Beta.idx[2],], index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai
        Betaitmp <- x$Betai
        etmp <- x$e
        sstmp <- x$ss

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,ix[i]]
            Betaitmp[,,i] <- x$Betai[,,ix[i]]
            etmp[,,i] <- x$e[,,ix[i]]
            sstmp$SS[,i] <- x$ss$SS[,ix[i]]
        }

        Q <- x$Q[ix,ix]
        sstmp$transitions <- x$ss$transitions[ix, ix]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp, ss=sstmp, e=etmp, Q=Q))
    }

    # Sort regimes based on values of the Sigma matrix element
    # selected.

    if(permute==FALSE & is.null(Sigma.idx)==FALSE){
        # Find the sort indices for the objects
        ix <- sort(x$Sigmai[Sigma.idx, Sigma.idx,], index.return=TRUE)$ix

        Sigmaitmp <- x$Sigmai
        Betaitmp <- x$Betai
        etmp <- x$e
        sstmp <- x$ss

        for(i in 1:h)
        {
            Sigmaitmp[,,i] <- x$Sigmai[,,ix[i]]
            Betaitmp[,,i] <- x$Betai[,,ix[i]]
            etmp[,,i] <- x$e[,,ix[i]]
            sstmp$SS[,i] <- x$ss$SS[,ix[i]]
        }

        Q <- x$Q[ix,ix]
        sstmp$transitions <- x$ss$transitions[ix, ix]
        return(list(Betai=Betaitmp, Sigmai=Sigmaitmp, ss=sstmp, e=etmp, Q=Q))
    }
}


####################################################################
# Full MCMC sampler for MSBVAR models with SZ prior
####################################################################

gibbs.msbvar <- function(x, N1=1000, N2=1000, permute=TRUE,
                         Beta.idx=NULL, Sigma.idx=NULL,
                         posterior.fit=FALSE)
{
    # Do some sanity / error checking on the inputs so people know
    # what they are getting back

    if(permute==TRUE & is.null(Beta.idx)==FALSE){
        cat("You have set permute=TRUE which means you requested a random permutation sample.\n Beta.idx argument will be ignored")
        }

    if(permute==TRUE & is.null(Sigma.idx)==FALSE){
        cat("You have set permute=TRUE which means you requested a random permutation sample.\n Sigma.idx argument will be ignored")
        }
    if(permute==FALSE & (is.null(Beta.idx)==TRUE &
        is.null(Sigma.idx)==TRUE)){
        cat("You have set permute=FALSE, but failed to provide an identification to label the regimes.\n Set Beta.idx or Sigma.idx != NULL")
    }

    # Initializations
    init.model <- x$init.model
    hreg <- x$hreg
    Q <- x$Q
    fp <- x$fp
    m <- x$m
    p <- x$p
    h <- x$h
    alpha.prior <- x$alpha.prior
    TT <- nrow(fp)

    # Initialize the Gibbs sampler parameters
    e <- hreg$e
    Sigmai <- hreg$Sigma
    Betai <- hreg$Bk

    # Burnin loop
    for(j in 1:N1)
    {
        # Smooth / draw the state-space
        ss <- SS.draw(e, TT, m, h, Sigmai, Q)

        # Draw Q
        Q <- MSBVAR.Q.draw(Q, ss$transitions, prior=alpha.prior)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, TT, ss$SS, init.model)

        # Draw the variance matrices for each regime
        Sigmai <- Sigma.draw(m, h, ss, hreg, Sigmai)

        # Sample the regression coefficients
        Betai <- Beta.draw(m, p, h, Sigmai, hreg, init.model, Betai)

        # Update the residuals
        e <- residual.update(m, h, init.model, Betai, e)

        # Now do the random permutation of the labels with an MH step
        # for the permutation or sort the regimes based on the
        # Beta.idx or Sigma.idx conditions.
        pf.obj <- PermuteFlip(x=list(Betai=Betai, Sigmai=Sigmai,
                                     ss=ss, e=e, Q=Q),
                              h, permute, Beta.idx, Sigma.idx)

        # Replace sampler objects with permuted / flipped ones
        Betai <- pf.obj$Betai
        Sigmai <- pf.obj$Sigmai
        ss <- pf.obj$ss
        e <- pf.obj$e
        Q <- pf.obj$Q

        # Print out some iteration information
        if(j%%1000==0) cat("Burn-in iteration : ", j, "\n")

    }

    # End of burnin loop!

    # Declare the storage for the return objects
    ss.storage <- vector("list", N2)
    transition.storage <- array(NA, c(h,h,N2))
    Beta.storage <- matrix(NA, N2, (m^2*p+m)*h)
    Sigma.storage <- matrix(NA, N2, m*(m+1)*0.5*h)
    Q.storage <- matrix(NA, N2, h^2)

    # Main loop after the burnin -- these are the sweeps that are kept
    # for the final posterior sample.
    # Note, we need to do the storage AFTER the permutation steps so
    # we get the labeling correct.
    for(j in 1:N2)
    {
        # Draw the state-space
        ss <- SS.draw(e, TT, m, h, Sigmai, Q)

        # Draw Q
        Q <- MSBVAR.Q.draw(Q, ss$transitions, prior=alpha.prior)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, TT, ss$SS, init.model)

        # Draw the variance matrices for each regime
        Sigmai <- Sigma.draw(m, h, ss, hreg, Sigmai)

        # Sample the regression coefficients
        Betai <- Beta.draw(m, p, h, Sigmai, hreg, init.model, Betai)

        # Update the residuals
        e <- residual.update(m, h, init.model, Betai, e)

        # Now do the random permutation of the labels with an MH step
        # for the permutation or sort the regimes based on the
        # Beta.idx or Sigma.idx conditions.
        pf.obj <- PermuteFlip(x=list(Betai=Betai, Sigmai=Sigmai,
                                     ss=ss, e=e, Q=Q),
                              h, permute, Beta.idx, Sigma.idx)

        # Replace sampler objects with permuted / flipped ones
        Betai <- pf.obj$Betai
        Sigmai <- pf.obj$Sigmai
        ss <- pf.obj$ss
        e <- pf.obj$e
        Q <- pf.obj$Q

        # Store things -- after flip / identification
        Sigma.storage[j,] <- as.vector(apply(Sigmai, 3, vech))
        Beta.storage[j,] <- as.vector(Betai)
        ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
        transition.storage[,,j] <- ss$transitions
        Q.storage[j,] <- as.vector(Q)

        # Print out some iteration information
        if(j%%1000==0) cat("Final iteration : ", j, "\n")
    }

    #################################################################
    # Compute the posterior fit if requested -- this is the logMDD
    # computation in the algo.
    #
    # This section should be later unloaded to posterior.fit with a
    # class.
    #################################################################

    pfit <- NA

    if(posterior.fit==TRUE & permute==TRUE)
    {
        cat("The options posterior.fit==TRUE and permute==TRUE do not make sense.\nInference is only meaningful for the identified regimes.\nPosterior fit measure has not been computed\n")
    }

    if(posterior.fit==TRUE & permute==FALSE)
    {

    # Find the posterior mode for each component of the model
    Beta.mean <- array(apply(Beta.storage, 2, mean), c(m*p + 1, m, h))
    Sigma.part <- matrix(apply(Sigma.storage, 2, mean), ncol=h)
    Sigma.mean <- array(sapply(1:h,
                               function(i){xpnd(Sigma.part[,i])}), c(m,m,h))

    Q.mean <- matrix(apply(Q.storage, 2, mean), h, h)

    ss.mean <- apply(matrix(unlist(lapply(ss.storage, as.integer)),
                            TT, N2), 1, sum)
    ss.mean <- matrix(ss.mean, TT, h-1)
    ss.mean <- cbind(ss.mean, rep(N2, TT) - rowSums(ss.mean))/N2

    # Compute the log priors
    B0 <- matrix(0, m*p + 1, m)
    diag(B0) <- 1
    BS <- kronecker(x$init.model$S0, solve(x$init.model$H0))
    log.B0 <- sum(sapply(1:h, function(i){dmvnorm(as.vector(Beta.mean[,,i]),
                                                  mean=as.vector(B0),
                                                  sigma=BS,
                                                  log=TRUE)}))

    log.S0 <- sum(sapply(1:h, function(i){ldwishart(solve(x$init.model$S0),
                                                    round(colSums(ss.mean))[i],
                                                    solve(Sigma.mean[,,i]))}))

    log.Q <- sum(log(ddirichlet(Q.mean, alpha.prior)))

    # Compute the log likelihood
    tmp <- hregime.reg2(h, m, p, TT, ss.mean, init.model)

    ll <- sum(sapply(1:h, function(i){dmvnorm(tmp$e[,,i],
                                              mean=rep(0,m),
                                              sigma=Sigma.mean[,,i],
                                              log=TRUE)}))

    # Find the pdf of the model parameters -- this is evaluating the
    # mode at each of the draws.

    pdf.Sigma <- pdf.Q <- rep(NA,N2)
    pdf.Beta <- matrix(NA, N2, h)

    hreg <- hregime.reg2(h, m, p, TT, ss.mean, init.model)

    # No need to permute here since we working with the identified
    # posterior mode.
    X <- init.model$X[(m+1):nrow(init.model$X),]

    for(i in 1:N2)
    {
        # Pr(Q | Beta, Sigma, y, S) -- resimulate the state-space,
        # redraw Q and then find the density at the mode.

        pdf.Q[i] <- prod(ddirichlet(matrix(Q.mean, h, h),
                                    transition.storage[,,i] + alpha.prior))

        # Pr(Sigma | Beta, Q, y, S)
        Sigmai <- matrix(Sigma.storage[i,], ncol=h)
        Sigmai <- array(sapply(1:h, function(j){xpnd(Sigmai[,j])}),
                        c(m,m,h))
        dftmp <- matrix(as.integer(ss.storage[[i]]), ncol=h-1)
        dftmp <- cbind(dftmp, 1 - rowSums(dftmp))
        df <- colSums(dftmp)
        pdf.Sigma[i] <- sum(sapply(1:h,
                                function(j){ldwishart(Sigma.mean[,,j],
                                                      df[j],
                                                      Sigmai[,,j])}))

        # Pr(Beta | Sigma, Q, y, S)

        # Compute the var-cov of the coefficients from the sample.
        # Need to get the X'X matrix for the state for the computation
        ss <- matrix(as.integer(ss.storage[[i]]), ncol=h-1)
        ss <- cbind(ss, rep(1, TT) - rowSums(ss))

        Betai <- matrix(Beta.storage[i,], ncol=h)

        for(j in 1:h)
        {
            XX <- crossprod(X, diag(ss[,j]))%*%X +
                  crossprod(init.model$X[1:(m+1),]) + init.model$H0

            BS <- kronecker(Sigmai[,,j], solve(XX))
            pdf.Beta[i,j] <- dmvnorm(as.vector(Beta.mean[,,j]),
                                     mean=Betai[,j],
                                     sigma=BS, log=TRUE)
        }

    }

    # log pdf components computation
    log.pdf.Q <- log(mean(pdf.Q))

    log.pdf.Sigma <- pdf.Sigma
    log.pdf.Sigma.max <- max(pdf.Sigma)
    qlog1 <- log.pdf.Sigma - log.pdf.Sigma.max
    log.pdf.Sigma <- log.pdf.Sigma.max - log(N2) + log(sum(exp(qlog1)))

    log.pdf.Beta <- apply(pdf.Beta, 1, sum)
    log.pdf.Beta.max <- max(log.pdf.Beta)
    qlog2 <- log.pdf.Beta - log.pdf.Beta.max
    log.pdf.Beta <- log.pdf.Beta.max - log(N2) + log(sum(exp(qlog2)))

    log.prior <- log.B0 + log.S0 + log.Q

    # Do the log marginal density calculation
    num <- ll + log.prior
    den <- log.pdf.Q + log.pdf.Sigma + log.pdf.Beta

    lmdd <- num - den

    # Store the results for the fit measures
    pfit <- list(log.llf=ll,
                 log.prior=log.prior,
                 log.marginal.data.density=lmdd,
                 log.pdf.Q=log.pdf.Q,
                 log.pdf.Sigma=log.pdf.Sigma,
                 log.pdf.Beta=log.pdf.Beta,
                 log.marginal.posterior=den)

}  # End of loop for the iterations for the posterior fit blocks

    # Now make an output object and provide classing
    class(ss.storage) <- c("SS")
    class(pfit) <- c("posterior.fit.MSBVAR")

    output <- list(Beta.sample=mcmc(Beta.storage),
                   Sigma.sample=mcmc(Sigma.storage),
                   Q.sample=mcmc(Q.storage),
                   transition.sample=transition.storage,
                   ss.sample=ss.storage,
                   pfit=pfit,
                   init.model=init.model,
                   h=h,
                   p=p,
                   m=m)

    class(output) <- c("MSBVAR")
    attr(output, "eqnames") <- attr(init.model, "eqnames")
    return(output)
}
