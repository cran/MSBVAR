# gibbs.msbvar.R Functions for Gibbs sampling an MSBVAR model based on
# a SZ prior.

# Patrick T. Brandt
# 20081113 : Initial version

# Draws of Q, the transition matrix
MSBVAR.Q.draw <- function(SStrans, prior)
{
    # compute the density based on the counts.
    alpha <- SStrans + prior
    Q <- t(apply(alpha, 1, rdirichlet, n=1))
    return(Q)
}

# MH sampler for Q
MSBVAR.Q.draw2 <- function(Q, SStrans, prior)
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
##     print(Q.old)
##     print(Q.new)
##     print(A)

    if( log(runif(1)) > A ) return(Q.old) else return(Q.new)
}

# vech function for efficiently sorting symmetric matrices (like Sigma)

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

# xpnd function for reconstructing symmetric matrices from vech

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

# State-space sampling

draw.SS1 <- function(u, TT, m, h, Sik, Q)
{
    u.vec <- matrix(0, (TT*m), h); sig.vec <- matrix(0,m*m,h)
    for(i in 1:h){
        u.vec[,i] <- as.vector(t(u[,,i]))
        sig.vec[,i] <- as.vector(Sik[,,i])
    }

    # Make a dummy object so we can reject degenerate samples
    ss <- list(SS=matrix(0, TT, h))

    # Call the sampler
    SQ <- steady.Q(Q)
##    cat("Drawing a state-space...")
     while(sum(colSums(ss$SS)<=1)>0)
     {

        ss <- .Call("SS.draw.cpp", u.vec, sig.vec, Q, SQ,
                    as.integer(TT), as.integer(h), as.integer(m))

    }
##    cat("completed a state-space.\n")
    return(ss)
}

draw.SS2 <- function(u, TT, m, h, Sik, Q)
{
    # Filter
    fp <- BHLK.filter(u, Sik, Q)
#    plot(ts(fp), plot.type="s", col=1:h)
    # Sample
    # Make a dummy object so we can reject degenerate samples
    ss <- list(SS=matrix(0, TT, h))
     while(sum(colSums(ss$SS)<=1)>0)
     {
        ss <- generate.states(fp, t(Q))
    }
    return(ss)

}

# Gibbs sampler for MSBVAR  models with SZ prior

gibbs.msbvar <- function(x, N1=1000, N2=1000, permute=TRUE, upper.idx=NULL, lower.idx=NULL)
{
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

    if(h>2) stop("Can only handle h=2 case right now!")

    # Initialize the Gibbs sampler
    e <- hreg$e
    Sigmai <- hreg$Sigma
    Betai <- hreg$Bk

    # Burnin loop
    for(j in 1:N1)
    {
        # smooth / draw
        ss <- draw.SS1(e, TT, m, h, Sigmai, Q)

        # Draw the Q
#        Q <- MSBVAR.Q.draw(ss$transitions, prior=alpha.prior)
        Q <- MSBVAR.Q.draw2(Q, ss$transitions, prior=alpha.prior)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, TT, ss$SS, init.model)

        # Draw the variances for the MSVAR from an inverse Wishart --
        # draws h matrices which are returned as an m^2 x h matrix.  These
        # then need to be unwound into the matrices we need.

        # Unscaled inverse Wishart draws
        tmp <- array(unlist(lapply(sapply(1:h,
                                          function(i) {rwishart(N=1,
                                                                df=sum(ss$SS[,i])+m,
                                                                Sigma=diag(m))},
                                          simplify=FALSE), solve)), c(m,m,h))

        df <- colSums(ss$SS)

        for (i in 1:h)
        {
            Sigmai[,,i] <- t(chol(hreg$Sigma[,,i]))%*%(tmp[,,i]*(df[i]))%*%chol(hreg$Sigma[,,i])
        }

       # Sample the regression coefficients and setup for the top of the
       # loop
        for(i in 1:h) {
            bcoefs.se <- t(chol(kronecker(Sigmai[,,i],
                                          solve(hreg$moment[[i]]$Sxx + init.model$hstar))))
            Betai[,,i] <- hreg$Bk[,,i] + matrix(bcoefs.se%*%matrix(rnorm(m^2*p+m), ncol=1), ncol=m)
            e[,,i] <- init.model$Y[(m+1):nrow(init.model$Y),] - init.model$X[(m+1):nrow(init.model$X),]%*%Betai[,,i]
        }

        # Now do the random permutation of the labels with an MH step
        # for the permutation
        if(permute==TRUE){
        rp <- rmultinom(1, h-1, rep(1/h, h))+1

        for(i in 1:h)
        {
            Sigmai[,,i] <- Sigmai[,,rp[i]]
            Betai[,,i] <- Betai[,,rp[i]]
            e[,,i] <- e[,,rp[i]]
            ss$SS[,i] <- ss$SS[,rp[i]]
        }

        Q <- Q[rp,rp]
        ss$transitions <- ss$transitions[rp, rp]
    }
        if(permute==FALSE & is.null(upper.idx)==FALSE &
           is.null(lower.idx)==FALSE)
         {
             tmp <- Betai
             if(Betai[upper.idx[1], upper.idx[2],1] >
                Betai[lower.idx[1], lower.idx[2],2])
            { next; } else {
                Betatmp <- Betai
                Sigmatmp <- Sigmai
                etmp <- e
                sstmp <- ss
                Qtmp <- Q

                # Do the flips
                Betai[,,1] <- Betatmp[,,2]
                Betai[,,2] <- Betatmp[,,1]

                Sigmai[,,1] <- Sigmatmp[,,2]
                Sigmai[,,2] <- Sigmatmp[,,1]

                e[,,1] <- etmp[,,2]
                e[,,2] <- etmp[,,1]

                ss$SS[,1] <- sstmp$SS[,2]
                ss$SS[,2] <- sstmp$SS[,1]

                ss$transitions[1,] <- rev(ss$transitions[2,])
                ss$transitions[2,] <- rev(ss$transitions[1,])

                Q[1,] <- rev(Qtmp[2,])
                Q[2,] <- rev(Qtmp[1,])

            }
         }
    }
    # End of burnin loop!

    # Declare the storage for the return objects
    ss.storage <- vector("list", N2)
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
        ss <- draw.SS1(e, TT, m, h, Sigmai, Q)

        # Draw the Q
#        Q <- MSBVAR.Q.draw(ss$transitions, prior=alpha.prior)
        Q <- MSBVAR.Q.draw2(Q, ss$transitions, prior=alpha.prior)

        # Update the regression step
        hreg <- hregime.reg2(h, m, p, TT, ss$SS, init.model)

        # Draw the variances for the MSVAR from an inverse Wishart --
        # draws h matrices which are returned as an m^2 x h matrix.  These
        # then need to be unwound into the matrices we need.

        # Unscaled inverse Wishart draws
        tmp <- array(unlist(lapply(sapply(1:h,
                                          function(i) {rwishart(N=1,
                                                                df=sum(ss$SS[,i])+m,
                                                                Sigma=diag(m))},
                                          simplify=FALSE), solve)), c(m,m,h))

        # Final draws from the Wisharts, scaled for the problem
        df <- colSums(ss$SS)
        for (i in 1:h)
        {
            Sigmai[,,i] <- t(chol(hreg$Sigmak[,,i]))%*%(tmp[,,i]*df[i])%*%chol(hreg$Sigmak[,,i])
        }

    # Sample the regression coefficients and setup for the top of the
    # loop
        for(i in 1:h) {
            bcoefs.se <- t(chol(kronecker(Sigmai[,,i],
                                          solve(hreg$moment[[i]]$Sxx + init.model$hstar))))
            Betai[,,i] <- hreg$Bk[,,i] + matrix(bcoefs.se%*%matrix(rnorm(m^2*p+m), ncol=1), ncol=m)
            e[,,i] <- init.model$Y[(m+1):nrow(init.model$Y),] - init.model$X[(m+1):nrow(init.model$X),]%*%Betai[,,i]

        }

        # Now do the random permutation of the labels with an MH step
        # for the permutation
        if(permute==TRUE){
        rp <- rmultinom(1, h-1, rep(1/h, h))+1

        for(i in 1:h)
        {
            Sigmai[,,i] <- Sigmai[,,rp[i]]
            Betai[,,i] <- Betai[,,rp[i]]
            e[,,i] <- e[,,rp[i]]
            ss$SS[,i] <- ss$SS[,rp[i]]
        }

        Q <- Q[rp,rp]
        ss$transitions <- ss$transitions[rp, rp]

        # Store things -- after flip
        Sigma.storage[j,] <- as.vector(apply(Sigmai, 3, vech))
        Beta.storage[j,] <- as.vector(Betai)
        ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
        Q.storage[j,] <- as.vector(Q)

        if(j%%1000==0) cat("Final iteration : ", j, "\n")

    }
        if(permute==FALSE & is.null(upper.idx)==FALSE &
           is.null(lower.idx)==FALSE)
         {
             tmp <- Betai
             if(Betai[upper.idx[1], upper.idx[2],1] >
                Betai[lower.idx[1], lower.idx[2],2])
             {
                # Store things -- no flip needed
                 Sigma.storage[j,] <- as.vector(apply(Sigmai, 3, vech))
                 Beta.storage[j,] <- as.vector(Betai)
                 ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
                 Q.storage[j,] <- as.vector(Q)

                 if(j%%1000==0) cat("Final iteration : ", j, "\n")
                 next;
           } else {
                Betatmp <- Betai
                Sigmatmp <- Sigmai
                etmp <- e
                sstmp <- ss
                Qtmp <- Q

                # Do the flips
                Betai[,,1] <- Betatmp[,,2]
                Betai[,,2] <- Betatmp[,,1]

                Sigmai[,,1] <- Sigmatmp[,,2]
                Sigmai[,,2] <- Sigmatmp[,,1]

                e[,,1] <- etmp[,,2]
                e[,,2] <- etmp[,,1]

                ss$SS[,1] <- sstmp$SS[,2]
                ss$SS[,2] <- sstmp$SS[,1]

                ss$transitions[1,] <- rev(ss$transitions[2,])
                ss$transitions[2,] <- rev(ss$transitions[1,])

                Q[1,] <- rev(Qtmp[2,])
                Q[2,] <- rev(Qtmp[1,])

                # Store things -- after flip
                Sigma.storage[j,] <- as.vector(apply(Sigmai, 3, vech))
                Beta.storage[j,] <- as.vector(Betai)
                ss.storage[[j]] <- as.bit.integer(as.integer(ss$SS[,1:(h-1)]))
                Q.storage[j,] <- as.vector(Q)

                 if(j%%1000==0) cat("Final iteration : ", j, "\n")
            }
         }
    }

   # Now make an output object
    class(ss.storage) <- c("SS")

    output <- list(Beta.sample=mcmc(Beta.storage),
                   Sigma.sample=mcmc(Sigma.storage),
                   Q.sample=mcmc(Q.storage),
                   ss.sample=ss.storage,
                   h=h)
    class(output) <- c("MSBVAR")
    return(output)
}

