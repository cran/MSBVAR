"posterior.fit" <- function(varobj, A0.posterior.obj=NULL)
{
    if(inherits(varobj, "VAR"))
    {

      stop("posterior.fit() not implemented for VAR class objects.\n")
##         output <- posterior.fit.VAR(varobj)
##         attr(output, "class") <- c("posterior.fit.VAR")
##         return(output)
    }
    if(inherits(varobj, "BVAR"))
    {
        output <- posterior.fit.BVAR(varobj)
        attr(output, "class") <- c("posterior.fit.BVAR")
        return(output)
    }
    if(inherits(varobj, "BSVAR"))
    {
        output <- posterior.fit.BSVAR(varobj, A0.posterior.obj)
        attr(output, "class") <- c("posterior.fit.BSVAR")
        return(output)
    }
}

"posterior.fit.BSVAR" <- function(varobj, A0.posterior.obj)
{  # Constants from the model object

    m <- ncol(varobj$A0.mode)
    gk <- varobj$F.posterior
    N2 <- A0.posterior.obj$N2

   # compute the pdf of the prior: p(A_0, A_+) and the exponential
   # term of pr(Y|A0, A+)

    yexpt <- 0.0              # initial value for the exp term of pr(Y|A0, A+)
    pa0ap <- matrix(0, m, 1)  # Vector of log priors.

    for (i in 1:m)
    {
      # Grab the parameters in each structural equation.
      a0k <- varobj$A0.mode[,i]
      apk <- varobj$F.posterior[,i]

      # covariances of the free parameters in each equation
      S0bar <- varobj$H0inv.tilde[[i]]
      Spbar <- varobj$Hpinv.tilde[,,i]

      # Next 2 lines are the free parameters in each equation
      bk <- t(varobj$Ui[[i]])%*%a0k
      gbark <- varobj$Pi.tilde[[i]]%*%bk

      # Compute the exponential term.
      yexpt <- yexpt + 0.5*(t(a0k)%*%varobj$YY%*%a0k -
                            2*t(apk) %*% varobj$XY %*% a0k +
                            crossprod(apk, varobj$XX)%*%apk)
      # log prior

      pa0ap[i,] <- -0.5*(ncol(matrix(varobj$Ui[[i]], m)) +
                         nrow(varobj$F.posterior))*log(2*pi) +
                             0.5*log(abs(sum(diag(S0bar)))) +
                                 0.5*log(abs(sum(diag(Spbar)))) -
                                     0.5*(t(bk)%*%S0bar%*%bk + t(gk[,i] - gbark)%*%Spbar%*%(gk[,i]-gbark))
  }


  # Find log p(Y | A0, A+) --- this is the log likelihood at
  # the restricted parameters, so it is conditional on the model.  We
  # can use this later for computing the MDD and log posterior for the
  # models.  We could compute this without the loop, but then we would
  # not have the cumulative version for the CBF calculations later!

#    log.llf <- matrix(0, varobj$df, 1)
    U <- qr.R(qr(t(varobj$A0.mode)))
    ada <- sum(log(abs(diag(U))))

    log.llf <- -0.5*m*varobj$df*log(2*pi) + varobj$df*ada + yexpt

##     for (i in 1:varobj$df)
##     {
##         log.llf.restricted[i,] <- -(0.5*m)*log(2*pi) + ada -
##             0.5*crossprod(varobj$structural.innovations[i,])
##     }
##     log.marginal.llf.restricted <- cumsum(log.marginal.llf.restricted)

    # Now find Pr(A+ | Y, A0)
    totparam <- length(varobj$F.posterior)
    ld.covar.param <- sum(sapply(1:m, function(i) {
      determinant(varobj$Hpinv.posterior[[i]])$modulus }))
    
    bk <- a2b(varobj$A0.mode, varobj$Ui)
    Apexpt <- 0
    n0cum <- c(0,cumsum(varobj$n0))

    for (i in 1:m)
    {
        bj <- bk[(n0cum[i]+1):(n0cum[(i+1)])]
        tmp <- varobj$F.posterior[,i] - varobj$P.posterior[[i]]%*%bj
        Apexpt <- Apexpt - 0.5%*%t(tmp)%*%varobj$Hpinv.posterior[[i]]%*%tmp
    }

    pa.ya0 <- -0.5*totparam*log(2*pi) + ld.covar.param + Apexpt

    ################################################################################
    # Now find the pdf for each equation: p(a0_k |Y, a0_{not k})
    ################################################################################
    # Set up the gibbs sampler
    setup <- gibbs.setup.bsvar(varobj)
    Tinv <- setup$Tinv
    UT <- setup$UT

    logpa0.yao <- matrix(0, m, 1)  # storage for log Pr(A0_k | Y A0_k*, k*\ne k)

    vlog <- matrix(0, N2, 1)       # storage across draws.
    A0gbs <- varobj$A0.mode
    Wout <- vector("list", m)
    bk <- a2b(varobj$A0.mode, varobj$Ui)
    b.free <- vector("list", m)
    for(i in 1:m)
    {
        b.free[[i]] <- bk[(n0cum[i]+1):(n0cum[(i+1)])]
    }

    # Build up the constant terms we need for each of the A0{k) calculations
    lgamma.term <- lgamma(0.5*(varobj$df+1))
    ldf <- log(varobj$df)
    df1 <- 0.5*(varobj$df+1)*log(2)
    constant1 <- 0.5*(varobj$df+varobj$n0)*ldf
    constant2 <- 0.5*(varobj$n0-1)*log(2*pi)
    constantterm <- constant1 - constant2 - df1 - lgamma.term

    logN2 <- log(N2)

    for (i in 1:(m-1))
    {
        # Loop over the N2 draws.....

          for(j in 1:N2)
          {
              Wout <- W.get(A0.posterior.obj$W.posterior, j)
              Vk <- solve(Tinv[[i]], Wout[[i]])
              gbeta <- solve(Vk, b.free[[i]])
              Vtr <- qr.R(qr(t(Vk)))

              vlog[j,1] <-
                  (constantterm[i] - log(abs(det(Vtr))) +
                   varobj$df*log(abs(gbeta[1])) -
                   0.5*varobj$df*crossprod(t(b.free[[i]]%*%solve(Vtr))))

          }

          # Compute the modified harmonic mean of the max.
          vlogmax <- max(vlog)
          qlog <- vlog-vlogmax
          vlogxhat <- vlogmax - logN2 + log(sum(exp(qlog)))
          logpa0.yao[i,1] <- vlogxhat

      }

    # For the last or m^th column

    drawi <- drawA0(A0gbs, UT, varobj$df, varobj$n0, Wout)
    A0gbs <- drawi$A0gbs
    Wout <- drawi$W
    Vk <- solve(Tinv[[m]], Wout[[m]])
    gbeta <- solve(Vk, b.free[[m]])
    Vtr <- qr.R(qr(t(Vk)))

    logpa0.yao[m,1] <-
        (constantterm[m] - log(abs(det(Vtr))) +
         varobj$df*log(abs(gbeta[1])) -
         0.5*varobj$df*crossprod(t(b.free[[m]]%*%solve(Vtr))))

    # Now compute the output objects

    log.prior <- sum(pa0ap)
    lPy <- log.prior + log.llf - sum(logpa0.yao) - pa.ya0

    output <- list(log.prior=log.prior,
                   log.llf=log.llf,
                   log.posterior.Aplus=pa.ya0,
                   log.marginal.data.density=lPy,
                   log.marginal.A0k=logpa0.yao)
    class(output) <- c("posterior.fit.BSVAR")
    return(output)
}

"posterior.fit.BVAR" <- function(varobj)
{
      # Assign variables
      capT <- varobj$pfit$capT
      m <- varobj$pfit$m
      ncoef <- varobj$pfit$ncoef
      num.exog <- varobj$pfit$num.exog
      nu <- varobj$pfit$nu
      H0 <- varobj$pfit$H0
      S0 <- varobj$pfit$S0
      Y <- varobj$pfit$Y
      X <- varobj$pfit$X
      hstar1 <- varobj$pfit$hstar1
      Sh <- varobj$pfit$Sh
      u <- varobj$pfit$u
      Bh <- varobj$pfit$Bh
      Sh1 <- varobj$pfit$Sh1

    # Compute the marginal posterior LL value
    # For the derivation of the integrand see Zellner 1971, Section 8.2

      scalefactor <- (sum(lgamma(nu + 1 - seq(1:m))) -
                      sum(lgamma(nu + capT + 1 - seq(1:m))))

    # Find some log-dets to make the final computation easier.
    # This does depend on the prior chosen, since some of these
    # matrices will be zero for flat prior model, so the ldet of
    # the S0 mtx will be zero.
    # This is done with if-else statements.

      M0 <- (diag(capT) + X%*%solve(H0)%*%t(X))
      B0 <- matrix(0,nrow=(ncoef+num.exog),ncol=m)
      diag(B0) <- 1
      Bdiff <- Y-X%*%B0

      ld.S0 <- determinant(S0, logarithm=T)
      ld.S0 <- ld.S0$sign * ld.S0$modulus

      ld.M0 <- determinant(M0, logarithm=T)
      ld.M0 <- ld.M0$sign * ld.M0$modulus

      ld.tmp <- determinant((S0 + t(Bdiff)%*%solve(M0)%*%Bdiff), logarithm=T)
      ld.tmp <- ld.tmp$sign * ld.tmp$modulus

      data.marg.llf <-  (- 0.5*capT*m*log(2*pi)
                         - m*0.5*ld.M0
                         + capT*0.5*ld.S0
                         - scalefactor
                         - 0.5*(nu+capT)*ld.tmp)

    # Now find the predictive posterior density
      M1 <- (diag(capT) + X%*%solve(hstar1)%*%t(X))

      ld.S1 <- determinant(Sh, logarithm=T)
      ld.S1 <- ld.S1$sign * ld.S1$modulus

      ld.M1 <- determinant(M1, logarithm=T)
      ld.M1 <- ld.M1$sign * ld.M1$modulus

      ld.tmp <- determinant((Sh + t(u)%*%solve(M1)%*%u), logarithm=T)
      ld.tmp <- ld.tmp$sign * ld.tmp$modulus

      data.marg.post <- (- 0.5*capT*m*log(2*pi)
                         - m*0.5*ld.M1
                         + capT*0.5*ld.S1
                         - scalefactor
                         - 0.5*(nu+capT)*ld.tmp)

    # Now compute the marginal llf and the posterior for the
    # coefficients
      Bdiff <- B0 - Bh
      ld.S1 <- determinant(Sh1, logarithm=T)
      ld.S1 <- ld.S1$sign * ld.S1$modulus
      wdof <- capT - ncoef - num.exog - m - 1

      scalefactor1 <- (wdof*m*0.5)*log(2) + 0.25*m*(m-1) + (sum(lgamma(wdof + 1 - seq(1:m))))
      scalefactor2 <- -0.5*(ncoef*m)*log(2*pi)
      coef.post <- (scalefactor1 + scalefactor2 -0.5*(nu + capT + m +1)*ld.S1
                    - 0.5*sum(diag(solve(Sh1)%*%Sh))
                    - 0.5*(ncoef+num.exog)*ld.S1
                    - 0.5*sum(diag(Sh1%*%t(Bdiff)%*%hstar1%*%Bdiff)))

      output <- list(data.marg.llf=data.marg.llf,
                     data.marg.post=data.marg.post,
                     coef.post=coef.post)
      class(output) <- c("posterior.fit.BVAR")
      return(output)
  }

## "posterior.fit.VAR" <- function(varobj)
## {
##       # Assign variables
##       capT <- varobj$pfit$capT
##       m <- varobj$m
##       ncoef <- varobj$pfit$ncoef
##       num.exog <- varobj$pfit$num.exog
##       nu <- 0
##       H0 <- 0*varobj$hstar
##       S0 <- 0*varobj$hstar
##       Y <- varobj$Y
##       X <- varobj$X
##       hstar1 <- varobj$hstar
##       Sh <- varobj$mean.S
##       u <- varobj$residuals
##       Bh <- varobj$Bhat
##       Sh1 <- varobj$pfit$Sh1

##       scalefactor <- (sum(lgamma(nu + 1 - seq(1:m))) -
##                       sum(lgamma(nu + capT + 1 - seq(1:m))))

##       M0 <- diag(capT)
##       B0 <- matrix(0,nrow=(ncoef+num.exog),ncol=m)
##       diag(B0) <- 1
##       Bdiff <- Y-X%*%B0

##       ld.tmp <- determinant((S0 + t(Bdiff)%*%solve(M0)%*%Bdiff), logarithm=T)
##       ld.tmp <- ld.tmp$sign * ld.tmp$modulus

##       data.marg.llf <-  (-0.5*capT*m*log(2*pi) - scalefactor -
##                          0.5*(nu+capT)*ld.tmp)

##       # Now find the predictive posterior density
##       M1 <- (diag(capT) + X%*%solve(hstar1)%*%t(X))

##       ld.S1 <- determinant(Sh, logarithm=T)
##       ld.S1 <- ld.S1$sign * ld.S1$modulus

##       ld.M1 <- determinant(M1, logarithm=T)
##       ld.M1 <- ld.M1$sign * ld.M1$modulus

##       ld.tmp <- determinant((Sh + t(u)%*%solve(M1)%*%u), logarithm=T)
##       ld.tmp <- ld.tmp$sign * ld.tmp$modulus

##       data.marg.post <- (-0.5*capT*m*log(2*pi) - m*0.5*ld.M1 +
##                          capT*0.5*ld.S1 - scalefactor  -
##                          0.5*(nu+capT)*ld.tmp)

##     # Now compute the marginal llf and the posterior for the
##     # coefficients
##       Bdiff <- B0 - Bh
##       ld.S1 <- determinant(Sh1, logarithm=T)
##       ld.S1 <- ld.S1$sign * ld.S1$modulus
##       wdof <- capT - ncoef - num.exog - m - 1

##       scalefactor1 <- (wdof*m*0.5)*log(2) + 0.25*m*(m-1) + (sum(lgamma(wdof + 1 - seq(1:m))))
##       scalefactor2 <- -0.5*(ncoef*m)*log(2*pi)
##       coef.post <- (scalefactor1 + scalefactor2 -0.5*(nu + capT + m +1)*ld.S1
##                     - 0.5*sum(diag(solve(Sh1)%*%Sh))
##                     - 0.5*(ncoef+num.exog)*ld.S1
##                     - 0.5*sum(diag(Sh1%*%t(Bdiff)%*%hstar1%*%Bdiff)))

##       output <- list(data.marg.llf=data.marg.llf,
##                      data.marg.post=data.marg.post,
##                      coef.post=coef.post)
##       class(output) <- c("posterior.fit.VAR")
##       return(output)
## }

"print.posterior.fit" <- function(x, ...)
{
##     if(inherits(x, "posterior.fit.VAR"))
##     { print.posterior.fit.VAR(x, ...) }

    if(inherits(x, "posterior.fit.BVAR"))
    { print.posterior.fit.BVAR(x, ...) }

    if(inherits(x, "posterior.fit.BSVAR"))
    { print.posterior.fit.BSVAR(x, ...) }
}

## "print.posterior.fit.VAR" <- function(x, ...)
## {
##     cat("Log marginal density, Pr(Y)         :  ", x$data.marg.llf, "\n")
##     cat("Predictive marginal density         :  ", x$data.marg.post, "\n")
##     cat("Coefficient marginal LLF/posterior  :  ", x$coef.post, "\n")
## }

"print.posterior.fit.BVAR" <- function(x, ...)
{
    cat("Log marginal density, Pr(Y)         :  ", x$data.marg.llf, "\n")
    cat("Predictive marginal density         :  ", x$data.marg.post, "\n")
    cat("Coefficient marginal LLF/posterior  :  ", x$coef.post, "\n")
}

"print.posterior.fit.BSVAR" <- function(x, ...)
{    cat("Log prior, Pr(A0, A+)            : ", x$log.prior, "\n")
     cat("Log LLF, Pr(Y | A0, A+)          : ", x$log.llf, "\n")
     cat("Log Pr(A+|Y, A0)                 : ", x$log.posterior.Aplus, "\n")
     cat("Log marginal density, Pr(Y)      : ", x$log.marginal.data.density, "\n")
     cat("Log marginal A0(k) | not A0(k)   : ", x$log.marginal.A0k, "\n")
}

