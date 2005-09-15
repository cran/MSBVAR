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
