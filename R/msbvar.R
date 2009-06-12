# msbvar() and related functions
# Patrick T. Brandt
# 20081113 : Initial version
#

# workhorse msbvar function with SZ prior.

msbvar <- function(y, z=NULL, p, h,
                   lambda0, lambda1, lambda3,
                   lambda4, lambda5, mu5, mu6, qm,
                   alpha.prior=100*diag(h) + matrix(2, h, h),
                   prior=0, max.iter=40)
{

    m <- ncol(y)

    # This should be part of a sanity.check.msbvar
    if(h==1)
    {
        stop("\n\n\t -- For MSBVAR models, h>1\n")
    }

    # Before the loop this is all initialization for the EM
    # implementation of the Bayesian model.

    # Now do a baseline, non-regime model using szbvar() since this
    # gives us all of the inputs we need for later
    init.model <- szbvar(Y=y, p, z=z, lambda0, lambda1, lambda3,
                         lambda4, lambda5, mu5, mu6, nu=ncol(y),
                         qm=4, prior=prior,
                         posterior.fit=FALSE)

    # Now set up an initial filter for the state-space
    # Get the residuals
    u <- init.model$residuals[(m+1):nrow(init.model$residuals),]
    u <- array(u, c(nrow(u), ncol(u), h))

    # Define an initial Q

    Q <- (1-(h*0.25/(h-1)))*diag(h) + matrix(0.25/(h-1), h, h)

    # Now initialize the filter / smoother / etc.  This is basically a
    # flat draw for the state probabilities

    TT <- nrow(u)
    fp <- rdirichlet(TT, rep(1,h))
    sp <- BHLK.smoother(fp, Q)
    fp1 <- sp[1,]
    fp[1,] <- fp1  # initial filter is the smoothed value!
    q.init <- Q[,1:(h-1)]

    # Loop for the EM algorithm: note we have already initiated an
    # E-step, so we can start with an M-step

    for(i in 1:max.iter)
    {
        # M-step
        # Compute Q
        Q <- count.transitions(round(sp)) + alpha.prior
        Q <- Q/rowSums(Q)

        # Regression step to get the regime-specific estimates
        hreg <- hregime.reg2(h, m, p, TT, sp, init.model)

        # E-step
        # Filter
        fp <- BHLK.filter(hreg$e, hreg$Sigma, Q)
        fp[1,] <- fp1  # reset 1st period to ML from the last smoother
        # Smoother
        sp <- BHLK.smoother(fp, Q)
        fp1 <- sp[1,]
    }

    # Do a final M-step hreg for the optimized values
    Q <- count.transitions(round(sp)) + alpha.prior
    Q <- Q/rowSums(Q)
    hreg <- hregime.reg2(h, m, p, TT, fp, init.model)

    # Set up the output
    output <- list(init.model=init.model,
                   hreg=hreg,
                   Q=Q,
                   fp=fp,
                   m=m, p=p, h=h,
                   alpha.prior=alpha.prior)

    class(output) <- c("MSVARsetup")
    return(output)
}


# This needs to be implemented in C // C++ code later
BHLK.smoother <- function(fp, Q)
{
    TT <- nrow(fp)
    h <- ncol(fp)
    p.smooth <- matrix(0, TT, h)

    p.smooth[TT,] <- fp[TT,]

    for(tt in (TT-1):1)
    {
        p.predict <- Q%*%fp[tt,]
        p.smooth[tt,] <- crossprod(Q, (p.smooth[tt+1,]/p.predict))*fp[tt,]
    }

    return(p.smooth)
}

hregime.setup <- function(h, m, p, TT)
{

    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- matrix(0, 1, h)
    e <- array(0, c(TT, m, h))
    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e))
}


# hregime.reg : computes moments for hregime regression for msbvar (as
# opposed to MSBSVAR model)  Note in this function we assume the
# regimes are observed (we round the regimes!)

hregime.reg <- function(h, m, p, TT, fp, Y, X, init.model)
{
    # Find obs. for each state
    ss <- matrix(0, nrow(fp), ncol(fp))
    ss[which(fp>1/h)] <- 1

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- matrix(0, 1, h)
    e <- array(0, c(TT, m, h))

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the matrix
        ytmp <- rbind(init.model$Y[1:(m+1),],
                      matrix(Y[which(ss[,i]==1),], ncol=m))
        xtmp <- rbind(init.model$X[1:(m+1),],
                        matrix(X[which(ss[,i]==1),], ncol=ncol(X)))

        Syy <- crossprod(ytmp) # *xi[j,i]
        Sxy <- crossprod(xtmp, ytmp) #*xi[j,i]
        Sxx <- crossprod(xtmp) #*xi[j,i]

        # Compute the regression coefficients

        hstar <- Sxx + init.model$H0
        Bk[,,i] <- solve(hstar,
                         (Sxy + init.model$H0[,1:m]))

        # Compute Sigma
        df[1,i] <- nrow(ytmp)-m-1
        Sigmak[,,i] <- (init.model$S0 + Syy + init.model$H0[1:m,1:m] -
                        t(Bk[,,i])%*%hstar%*%Bk[,,i])/(df[1,i])

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        # Save the moments
        tmp[[i]] <- list(Syy=Syy, Sxy=Sxy, Sxx=Sxx, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}



# Regime weighted regression function for msbvar() -- this does the
# M-step for the regimes-dependent MSBVARs.  This version lets
# everything vary (intercepts, AR, Sigma).  Need to later implement
# versins ala Krolzig that are variations on this.

hregime.reg2 <- function(h, m, p, TT, fp, init.model)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- colSums(fp)
    e <- array(0, c(TT, m, h))
    Y <- init.model$Y[(m+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1):nrow(init.model$X),]

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices
        Sxy <- crossprod(X, diag(fp[,i]))%*%Y + crossprod(init.model$X[1:(m+1),], init.model$Y[1:(m+1),])
        Sxx <- crossprod(X, diag(fp[,i]))%*%X + crossprod(init.model$X[1:(m+1),])

        # Compute the regression coefficients
        hstar <- Sxx + init.model$H0
        Bk[,,i] <- solve(hstar,
                         (Sxy + init.model$H0[,1:m]))

        # Compute residuals and Sigma (based on Krolzig)

##         Sigmak[,,i] <- (init.model$S0 + Syy + init.model$H0[1:m,1:m] -
##                         t(Bk[,,i])%*%hstar%*%Bk[,,i])/(df[i])

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        Sigmak[,,i] <- (init.model$S0 + crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]

        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx) #, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}

hregime.reg3 <- function(h, m, p, TT, fp, init.model, hold)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- hold$Bk
    Sigmak <- hold$Sigmak
    e <- hold$e

    df <- hold$df
    Y <- init.model$Y[(m+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1):nrow(init.model$X),]

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices
        Sxy <- crossprod(X, diag(fp[,i]))%*%Y + crossprod(init.model$X[1:(m+1),], init.model$Y[1:(m+1),])
        Sxx <- crossprod(X, diag(fp[,i]))%*%X + crossprod(init.model$X[1:(m+1),])

        # Compute the regression coefficients
        hstar <- Sxx + init.model$H0
##         Bk[,,i] <- solve(hstar,
##                          (Sxy + init.model$H0[,1:m]))

        Bk[,,i] <- solve(hstar, (Sxy + init.model$H0[,1:m]))

        # Compute residuals and Sigma (based on Krolzig)

##         Sigmak[,,i] <- (init.model$S0 + Syy + init.model$H0[1:m,1:m] -
##                         t(Bk[,,i])%*%hstar%*%Bk[,,i])/(df[i])

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        Sigmak[,,i] <- (init.model$S0 + crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]

#        Sigmak[[i]] <- (init.model$S0 + crossprod(e[[i]],diag(fp[,i]))%*%e[[i]])/df[i]
        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx) #, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}

# cut this later!
count.transitions <- function(s)
  { M <- ncol(s)
    TT <- nrow(s)
    s <- crossprod(t(s), as.matrix(seq(1:M)))
    sw <- matrix(0, M, M)
    for (t in 2:TT)
      { st1 <- s[t-1]
        st <- s[t]
        sw[st1,st] <- sw[st1, st] + 1
      }
    return(sw)
  }


## # MSBVAR Q.posterior objective function for optimization

## MSVAR.Q.posterior <- function(q.init, alpha.prior, h, fp)
## {
##     Q.full <- matrix(cbind(matrix((q.init), h, h-1),
##                            matrix(1, h, 1) - rowSums(matrix(q.init, h, h-1))), h, h)

## #    fp <- BHLK.filter(hreg$e, hreg$Sigma, Q.full)
## #    fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize
##     ct <- count.transitions(round(fp))
##     alpha <- ct + alpha.prior - matrix(1, h, h)
##     fn <- sum(log(ddirichlet(Q.full, alpha)))
##     return(fn)

## }

## # MSVAR.Q.posterior gradient function

## MSVAR.Q.grad <- function(q.init, alpha.prior, h, fp)
## {
##     Q.full <- matrix(cbind(matrix((q.init), h, h-1),
##                            matrix(1, h, 1) - rowSums(matrix(q.init, h, h-1))), h, h)

## #    fp <- BHLK.filter(hreg$e, hreg$Sigma, Q.full)
## #    fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize
##     ct <- count.transitions(round(fp))
##     alpha <- ct + alpha.prior - matrix(1, h, h)
##     return(q.init/as.vector(alpha[,1:(h-1)]))

## }

