# hidden.R
# Patrick T. Brandt
# brandt@unt.edu

# This file contains utility functions for VAR models.  These
# functions are not really meant to be called by users.  In general,
# they reformulate and reshape the VAR parameters so that they are
# easier to work with.

# Changes in these functions change their applications in all related
# files in the package.

# Compute VAR forecasts using the coefficients that have been parsed
# from the VAR object.
"coef.forecast.var" <-
function(y, intercept, ar.coefs, exog.coefs, m, p,
                       capT, nsteps, A0=t(chol(var.obj$mean.S)),
                       shocks=matrix(0,nrow=nsteps,ncol=m),
                       exog.fut=matrix(0,nrow=nsteps, ncol=nrow(exog.coefs)))
{
  yhat<-rbind(y,matrix(0,ncol=m,nrow=nsteps))

   # Compute the deterministic part of the forecasts (less the intercept!)
   if(is.na(sum(exog.coefs))==F)
     {
       deterministic.var <- as.matrix(exog.fut) %*% exog.coefs
     }
   else
     { deterministic.var <- matrix(0,nrow=nsteps,ncol=m)
     }

   # Now loop over the forecast horizon
   for(h in 1:nsteps)
     {  yhat[capT + h, ] <- (yhat[capT + h - 1,] %*% ar.coefs[,,1] +
                             intercept + deterministic.var[h,] + (shocks[h,]%*%A0))
       if (p>1) {for(i in 2:p)
       { yhat[capT + h, ] <- (yhat[capT + h, ] +
                              (yhat[capT + h - i, ] %*% ar.coefs[,,i]))

       }}
     }
   return(ts(yhat))
 }

# Compute impulse responses using vec(ar.coefs | exog coefs).  This is
# a utility function, since it speeds computation of the VAR IRFs in
# Monte Carlo simulations where the vec(coefs) is drawn.
"irf.var.from.beta" <-
function(A0,bvec,nsteps)
  { m <- ncol(A0)
    p <- (length(bvec))/m^2;
    bmtx <- matrix(bvec,ncol=m)
    ar.coefs<-t(bmtx)     # extract the ar coefficients
    dim(ar.coefs)<-c(m,m,p)              # push ar coefs into M x M x P array
    ar.coefs<-aperm(ar.coefs,c(2,1,3))   # reorder array so columns
                                         # are for eqn

    impulses <- irf.var(list(ar.coefs=ar.coefs),nsteps,A0=A0)$mhat
    impulses <- aperm(impulses, c(3,1,2)) # flips around the responses
                                          # to stack them as series
                                          # for each response, so the
                                          # first nstep elements are
                                          # the responses to the first
                                          # shock, the second are the
                                          # responses of the first
                                          # variable to the next
                                          # shock, etc.

    dim(impulses)<-c((m^2)*nsteps,1)

    return(impulses)
  }

# posterior fit measures for szbvar
"posterior.fit.szbvar" <-
function(capT,m,ncoef,num.exog,nu,H0,S0,Y,X,hstar1,Sh,u, Bh,Sh1)
  {
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


    return(list(data.marg.llf=data.marg.llf,
                data.marg.post=data.marg.post,
                coef.post=coef.post))
  }


# Marginal Log-Posterior for A0 -- returns a single value.  Used for
# optimization of the peak of the posterior pdf for B-SVAR model.

# Inputs are part of the moments generated by STRUCTURAL.SZBVAR.
# b = free parameter vector in A0
# Ui = U matrices that map from a to b
# df = degrees of freedom : T - p + ndum
# H0inv.posterior = m-array of the covariances for the columns of the
# A0 matrix.

"A0.llf" <- function(b, Ui, df, H0inv.posterior)
{     m <- length(H0inv.posterior)
      n0 <- sapply(1:m, function(i) {ncol(as.matrix(Ui[[i]]))})
      n0cum <- c(0,cumsum(n0))
      llf <- 0

      A0 <- matrix(0, m, m)
      for (i in 1:m)
      {
          bj <- b[(n0cum[i]+1):(n0cum[(i+1)])]
          A0[, i] <- as.matrix(Ui[[i]])%*%bj
          llf <- llf + 0.5*(t(bj)%*%H0inv.posterior[[i]]%*%bj)
      }

      tmp <- lu(Matrix(t(A0), sparse=FALSE))
      L <- matrix(expand(tmp)$L@x, m, m)
      U <- matrix(expand(tmp)$U@x, m, m)
      llf <- llf - df*sum(log(abs(diag(U))))
      return(llf)
}

# TRANSLATION FUNCTIONS FOR SQUEEZED PARAMETERS IN structural B-SVAR
# models:

# Some utility functions for b2a and a2b.  These convert from the
# squeezed parameter vector b_i to a_i using b_i = U_i' a_i as n
# footnote 11 of W&Z 2003.

"b2a" <- function(b, Ui)
{ m <- nrow(Ui[[1]])
  n0 <- sapply(1:m, function(i) {ncol(as.matrix(Ui[[i]]))})
  n0cum <- c(0,cumsum(n0))
  A0 <- matrix(0,m,m)

  for(i in 1:m)
    {
      A0[,i] <- as.matrix(Ui[[i]])%*%b[(n0cum[i]+1):(n0cum[(i+1)])]
    }

  return(A0)
}

# A0 = full A0 matrix....
"a2b" <- function(A0, Ui)
  { m <- nrow(A0)
    n0 <- sapply(1:m, function(i) {ncol(as.matrix(Ui[[i]]))})
    n0cum <- c(0,cumsum(n0))
    b <- matrix(0, n0cum[length(n0cum)], 1)

    for(i in 1:m)
      {
        b[(n0cum[i]+1):(n0cum[(i+1)])] <- t(as.matrix(Ui[[i]]))%*%A0[,i]
      }
    return(b)
  }


# Set up function for the Gibbs sampler for A0 for structural.szbvar
# and gibbs.A0.  Input is a fitted structural.szbvar() object.

"gibbs.setup.bsvar" <- function(szbvar.obj)
  {  # Compute some objects we need for inputs for the Gibbs sampler
     # for the B-SVAR models
    m <- ncol(szbvar.obj$A0.mode)
     # T_i* = decomposition of the qi x qi covariance of the squeezed A0
     #        parameters for each equation

    Tinv <- sapply(1:m, function(i)
                   {chol(szbvar.obj$H0inv.posterior[[i]])/szbvar.obj$df},
                   simplify=F)

    # UT = m x qi matriX in equation 14 of WZ JEDC --- eqn 14.
    UT <- sapply(1:m, function(i)
                 {szbvar.obj$Ui[[i]]%*%solve(Tinv[[i]])},
                 simplify=F)

   ##  # inverse cholesky of Hpinv.posterior
##     VHphalf <- sapply(1:m, function(i)
##                       {solve(chol(H0inv.posterior[[i]]))}, simplify=F)
##     # Squeezed posterior of the parameters
##     PU <- sapply(1:m, function(i) { P.posterior[[i]]%*%t(Ui[[i]]) }, simplify=F )
##     VPU <- PU
##
    return(list(UT=UT, Tinv=Tinv))
  }


# ONE STEP Function to draw 1 A0 from the posterior of a structural
# BVAR model.  Can be used for the over and just identified cases.
# Inputs all come from setup function above or the model object.
"drawA0" <- function(A0gbs, UT, df, n0, Wout)
  { # Set up constants and place holder matrices
    m <- dim(A0gbs)[1]

    for(i in 1:m)
      {
        w <- matrix(0, m, 1)

        X <- A0gbs
        X[,i] <- 0
        tmp <- lu(Matrix(t(X), sparse=FALSE))

        L <- matrix(expand(tmp)$L@x, m, m)
        U <- matrix(expand(tmp)$U@x, m, m)
        d.U <- abs(diag(U))
        sing.values <- d.U/max(d.U)
        jIx0 <- min(which(sing.values<2^-30))
        w[jIx0] <- 1
        if(jIx0>1 && jIx0<(m+1))
          {
            jA <- U[1:(jIx0-1), 1:(jIx0-1)]
            jb <- U[1:(jIx0-1), jIx0]
            jy <- qr.solve(-jA, jb)
            w[1:(jIx0-1)] <- jy
          }

        # now find the orthonormal basis for the w_i... w_qi
        w0 <- crossprod(UT[[i]], w)
        w1 <- w0/sqrt(sum(w0^2))

     # This is a bit of a "cheat", but it lets us trap occurences
     # where the solution is degenerate or the QR solution is fragile.  The
     # w converge quickly (within the first 100 iterations, so reusing the
     # previous ones in this error trap is probably OK.  If nothing
     # else, it increases the numerical stability of the algorithm

        if(is.na(sum(w1))==TRUE)
        {    W <- Wout[[i]]
          } else {
              W <- qr.Q(qr(w1, tol=1e-12), complete=T)
          }

        if(is.numeric(W)==FALSE) W <- Wout[[i]]

        gkb <- matrix(0, n0[i], 1)
        jstd <- sqrt(1/df)
        if(n0[i]>1)
          {
            gkb[2:n0[i],1] <- jstd*rnorm(n0[i]-1)
          }
        jr <- jstd*rnorm(df+1)
        gkb[1,1] <- ifelse(runif(1)<0.5, sqrt(crossprod(jr)),
                           -1*sqrt(crossprod(jr)))
        A0gbs[,i] <- UT[[i]]%*%(W%*%gkb)
        Wout[[i]] <- W
      }
    return(list(A0gbs=A0gbs, W=Wout))

  }

