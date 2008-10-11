"forecast" <- function(varobj, nsteps, A0,
                       shocks=matrix(0,nrow=nsteps,ncol=dim(varobj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps,ncol=nrow(varobj$exog.coefs)))
{
    if(inherits(varobj,"VAR")){
        return(forecast.VAR(varobj, nsteps, A0=t(chol(varobj$mean.S)),
                            shocks=shocks, exog.fut=exog.fut))
    }

    if(inherits(varobj, "BVAR")){
        return(forecast.VAR(varobj, nsteps, A0=t(chol(varobj$mean.S)),
                            shocks=shocks, exog.fut=exog.fut))
    }

    if(inherits(varobj, "BSVAR")){
        return(forecast.VAR(varobj, nsteps, A0=solve(varobj$A0.mode),
                       shocks=shocks, exog.fut=exog.fut))
    }
}

"forecast.VAR" <-
function(varobj, nsteps, A0=t(chol(varobj$mean.S)),
                       shocks=matrix(0,nrow=nsteps,ncol=dim(varobj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps,ncol=nrow(varobj$exog.coefs)))
{
   # Set up the initial parameters for the VAR forecast function from
   #  VAR object
  y <- varobj$y
  intercept <- varobj$intercept
  ar.coefs <- varobj$ar.coefs
  exog.coefs <- varobj$exog.coefs
  m<-dim(ar.coefs)[1]
  p<-dim(ar.coefs)[3]
  capT<-nrow(y)
  yhat<-rbind(y,matrix(0,ncol=m,nrow=nsteps))

   # Compute the deterministic part of the forecasts (less the intercept!)
   if(is.na(sum(varobj$exog.coefs))==F)
     {
       deterministic.VAR <- as.matrix(exog.fut) %*% exog.coefs
     }
   else
     { deterministic.VAR <- matrix(0,nrow=nsteps,ncol=m)
     }

   # Now loop over the forecast horizon
   for(h in 1:nsteps)
     {  yhat[capT + h, ] <- (yhat[capT + h - 1,] %*% ar.coefs[,,1] +
                             intercept + deterministic.VAR[h,] + (shocks[h,]%*%A0))
       if (p>1) {for(i in 2:p)
       { yhat[capT + h, ] <- (yhat[capT + h, ] +
                              (yhat[capT + h - i, ] %*% ar.coefs[,,i]))

       }}
     }
  output <- ts(yhat)
  attr(output, "class") <- c("forecast.VAR", "mts", "ts")
  return(output)
}

"forecast.BVAR" <- function(varobj, nsteps, A0=t(chol(varobj$mean.S)),
                       shocks=matrix(0,nrow=nsteps,ncol=dim(varobj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps,ncol=nrow(varobj$exog.coefs)))
{
    output <- forecast.VAR(varobj, nsteps, A0, shocks, exog.fut)
    attr(output, "class") <- c("forecast.BVAR", "mts", "ts")
    return(output)
}

"forecast.BSVAR" <- function(varobj, nsteps, A0=solve(varobj$A0.mode),
                       shocks=matrix(0,nrow=nsteps,ncol=dim(varobj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps,ncol=nrow(varobj$exog.coefs)))
{
    output <- forecast.VAR(varobj, nsteps, A0, shocks, exog.fut)
    attr(output, "class") <- c("forecast.BSVAR", "mts", "ts")
    return(output)
}

"uc.forecast" <- function(varobj, nsteps, burnin, gibbs, exog=NULL)
{
    if(inherits(varobj, "VAR"))
    {
        output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BVAR"))
    {
        output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BSVAR"))
    {
      stop("Not yet implemented for BSVAR models!\n")
##       output <- uc.forecast.VAR(varobj, nsteps, burnin, gibbs,exog)
##       attr(output, "class") <- c("uc.forecast.VAR", "mts", "ts")
##       return(output)
    }
}

"uc.forecast.VAR" <- function(varobj, nsteps, burnin, gibbs, exog=NULL)
  { # Extract all the elements from the VAR object
    y <- varobj$y
    ar.coefs <- varobj$ar.coefs
    intercept <- varobj$intercept
    A0 <- t(chol(varobj$mean.S))
    X <- varobj$X         # rhs variables for the model
    Y <- varobj$Y         # lhs variables for the model
    H0 <- varobj$H0
    S0 <- varobj$S0
    mu <- varobj$hyperp
    exog.coefs <- varobj$exog.coefs
    z <- varobj$z
    lambda0 <- varobj$prior[1]
    lambda1 <- varobj$prior[2]
    lambda3 <- varobj$prior[3]
    lambda4 <- varobj$prior[4]
    lambda5 <- varobj$prior[5]
    mu5 <- varobj$prior[6]
    mu6 <- varobj$prior[7]
    nu <- varobj$prior[8]
    prior <- varobj$prior.type
    num.exog <- varobj$num.exog
    qm <- varobj$qm
    ncoef <- nrow(varobj$Bhat)

  # Get some constants we are going to need
    starttime<-date()           # Starting time for simulation

    p<-dim(ar.coefs)[3]         # Capture the number of lags
                                # from input ar coefficients

    m<-ncol(y);                 # Number of endogenous
                                # variables in the VAR
    k<-m*nsteps;              # k = mh, the maximal number
                              # of forecasts

    capT<-nrow(y)               # Number of observations we
                                # are going to use.

  # Make arrays to hold the Gibbs sampler results
  yforc<-array(0,c(gibbs,nsteps,m))

  # Do the Gibbs draws.....
  for(i in 1:(burnin+gibbs))
    { # Step (a): Compute a draw of the conditional forecasts
      # COMPUTE INNOVATIONS: These are the structural innovations

      # First draw the innovations
      epsilon.i <- matrix(rnorm(nsteps*m),nrow=nsteps,ncol=m)

      # Then construct a forecast using the innovations
      ytmp <- forecast.VAR(varobj, nsteps, A0=A0, shocks=epsilon.i)

      # Store draws that are past the burnin in the array
      if(i>burnin)
        { for(j in 1:m)
            { yforc[(i-burnin),(1:nsteps),j]<-ytmp[((capT+1):(capT+nsteps)),j] }
        }


      # Step (b): Compute the mode of the posterior for the
      # forecast distribution.  This is the "extended"
      # dataset that includes the i'th Gibbs sample forecast

      # Set up the updated Y Matrix
      # this is just "ytmp" from above

      Y.update <- ytmp[(capT-p+1):nrow(ytmp),]

      # Set up the updated X -- this is hard because we need to get
      # the RHS lags correct.  We do this by padding the existing Y
      # and then building the lags.  This reuses the code for the lag
      # construction in the szbvar code

      # Now build rhs -- start with an empty matrix
      X.update <- matrix(0, nsteps, m*p+1)
        # Add in the constant
      X.update[, (m*p+1)] <- matrix(1, nsteps, 1)
        # Put in the lagged y's

      # Note that constant is put last here when the lags are made
      for(j in 1:p)
        {
          X.update[1:nsteps,(m*(j-1)+1):(m*j)] <- matrix(Y.update[(p+1-j):(nsteps+p-j),],ncol=m)
        }

      # Put in exogenous coefficients if there are any.
      if(is.null(exog)==F)
        {
          X.update<-cbind(X.update,exog);
        }

      # Now, stack the original Y data and the augmented data.
      Y.update <- rbind(Y, ytmp[(capT+1):nrow(ytmp),])

      # Set up crossproducts and inverses we need
      X.update <- rbind(X, X.update)
      XX.update <- crossprod(X.update)     # Cross product of RHS variables
      hstar.update <- H0 + XX.update       # Prior + Cross product

      # Updated Regression estimates, Beta | Sigma
      B.update<-solve((hstar.update),(crossprod(X.update,Y.update) + H0[,1:m]))

      # Posterior mean of Sigma | Beta
      S.update <- (S0 + crossprod(Y.update)
                   + H0[1:m,1:m] - t(B.update)%*%(hstar.update)%*%(B.update))/(capT+nsteps+nu-m-1)

      # Posterior variance of B: VBh = diag(Sh.*.(inv((H0 + x'x)))
      hstarinv <- solve(hstar.update)
      vcv.Bh <- kronecker(S.update,hstarinv)

      # Draw from the conditional posterior pdfs of the parameters

      # This is only valid for just-identified models.
      df <- capT - m*p - m - 1 + nsteps
      wisharts <- rwishart(1, df, diag(m))

      # Generate the draws from the Wishart and the Beta
      # Wishart draw
      Sigmat <- (chol(S.update))
      Sigma.Draw <- t(Sigmat)%*%(df*solve(matrix(wisharts,m,m)))%*%Sigmat
      sqrtwish <- t(chol(Sigma.Draw))
      # Covariance of beta
      bcoefs.covar <- t(chol(vcv.Bh))

      # Draw beta|Sigma
      aplus <- matrix(B.update, ncol=1) + bcoefs.covar%*%matrix(rnorm(nrow(bcoefs.covar)), ncol=1)
      aplus <- matrix(aplus, ncol=m)

      aplus.coefs<-t(aplus[1:(m*p),]);        # extract the ar coefficients
      dim(aplus.coefs)<-c(m,m,p)                    # push ar coefs into M x M x P array
      aplus.coefs<-aperm(aplus.coefs,c(2,1,3))      # reorder array


      intercept <- aplus[(m*p+1),]       # get drawn intercept....
      ar.coefs<-aplus.coefs            # AR coefs
      A0 <- sqrtwish

      if(num.exog!=0)
        {
          exog.coefs <- aplus[(m*p+2):nrow(aplus),]
        }


      # Print some intermediate results to capture progress....
      # and tell us that things are still running
      if (i%%500==0)
        { cat("Gibbs Iteration = ", i, "     \n");
          if(i<=burnin)
            { cat("(Still a burn-in draw.)\n");
            }

        }
      # Back to the top of the Gibbs loop....
    }
  endtime<-date()
  # Print time stamp so we know how long everything took.
  cat("Start time : ", starttime, "\n");
  cat("End time   : ", endtime, "\n");
    # Returns a list object
    output <- list(forecast=yforc)
    attr(output, "class") <- c("forecast.VAR")
    return(output)
}

"hc.forecast" <- function(varobj, yconst, nsteps, burnin, gibbs, exog=NULL)
{
    if(inherits(varobj, "VAR"))
    {
        output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
                                  gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BVAR"))
    {
        output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
                                  gibbs, exog)
        attr(output, "class") <- c("forecast.VAR")
        return(output)
    }

    if(inherits(varobj, "BSVAR"))
    {
      stop("Not yet implemented for B-SVAR models!\n")
##         output <- hc.forecast.VAR(varobj, yconst, nsteps, burnin,
##                                   gibbs, exog)
##         attr(output, "class") <- c("hc.forecast.VAR", "mts", "ts")
##         return(output)
    }
}

"hc.forecast.VAR" <-
function(varobj, yconst, nsteps, burnin, gibbs, exog=NULL)
{
    # Extract all the elements from the VAR object that we will need
    y <- varobj$y
    ar.coefs <- varobj$ar.coefs
    intercept <- varobj$intercept
    exog.coefs <- varobj$exog.coefs
    A0 <- t(chol(varobj$mean.S))
    mu <- varobj$hyperp
    prior <- varobj$prior
    ncoef <- nrow(varobj$Bhat)
    X <- varobj$X         # rhs variables for the model
    Y <- varobj$Y         # lhs variables for the model
    H0 <- varobj$H0       # precision for the var coefs
    S0 <- varobj$S0       # precision for the Sigma
    nu <- varobj$prior[8] # df

    # Get some constants we are going to need from the inputs
    starttime<-date()           # Starting time for simulation
                                # of forecasts

    q<-nrow(as.matrix(yconst));   # Number of restrictions

    p<-dim(ar.coefs)[3]         # Capture the number of lags
                                # from input ar coefficients

    m<-ncol(y);                 # Number of endogenous
                                # variables in the VAR

    capT<-nrow(y)               # Number of observations we
                                # are going to use.

    k<-m*nsteps;                # k = mh, the maximal number

    q <- nrow(yconst)             # Number of constraints

    # Make arrays to hold the Gibbs sampler results
    yforc<-array(0,c(gibbs,nsteps,m))


    # Do the Gibbs draws.....
    for(i in 1:(burnin+gibbs))
    {
      # Step (a): Compute a draw of the conditional forecasts
      # COMPUTE INNOVATIONS: These are the structural innovations
      # Solve the constraint equation for the updated forecast errors

      # Generate the forecasts without shocks
      ytmp<-as.matrix(coef.forecast.VAR(y, intercept, ar.coefs, exog.coefs, m, p,
                                     capT, nsteps, A0=(A0)))

      # Get the impulse responses that correspond to the forecasted data
      M <- irf.VAR(varobj, nsteps, A0)$mhat

      # Construct the draw of the orthogonalized innovations that
      # satisfy the hard condition.

      # These are the q constrained innovations
      r <- (yconst - ytmp[(capT+1):(capT+nsteps),])
      r<-matrix(r[1:nsteps,1],ncol=1)

      # Build the matrix of the impulses that define the constraint
      R <- matrix(0, k, q)

      # Put the g'th column of the impulse into the constraint matrix,
      # such that R * epsilon = r
      for (g in 1:q)
      {
        if(g==1)
          { R[1:length(M[,1,1]), 1] <- M[,1,1] }
        else
          {
            R[,g] <- c(M[,1,g], R[,g-1])[1:k]
          }
      }

      # Solve the minimization problem for the mean and variance of
      # the constrained innovations.

      RRinv<-solve(crossprod(R))
      mean.epsilon <- R%*%RRinv%*%r;
      var.epsilon <- diag(1,nrow=k) - (R%*%RRinv%*%t(R));

      # Draw from the singular MVN pdf of the constrained innovations.

      epsilon.i <- matrix(rmultnorm(1, mean.epsilon, var.epsilon),
                          nrow=nsteps, ncol=m, byrow=T)

      # Add the innovations to the forecasts
      ytmp[(capT+1):(capT +nsteps),] <- ytmp[(capT+1):(capT +nsteps),]+epsilon.i%*%A0

      # Store forecasts that are past the burnin point
      if(i>burnin)
        { for(j in 1:m)
            { yforc[(i-burnin),(1:nsteps),j]<-ytmp[(capT+1):(capT+nsteps),j] }
        }


      # Step (b): Compute the mode of the posterior for the
      # conditional forecast distribution.  This is the "extended"
      # dataset that includes the i'th Gibbs sample forecast


      # Build the augmented LHS and RHS matrices
      # 1) Get the nsteps+p observations we need to build the lagged
      # endogenous variables for the augmented system.

      Y.update <- ytmp[(capT-p+1):nrow(ytmp),]

      # 2) Build the updated X -- this is hard because we need to get
      # the RHS lags correct.  We do this by padding the existing Y
      # and then building the lags.  This reuses the code for the lag
      # construction in the szbvar code

      X.update <- matrix(0, nsteps, ncoef)
      X.update[,ncoef] <- matrix(1, nsteps, 1)

      # Note that constant is put last here when the lags are made
      for(j in 1:p)
        {
          X.update[1:nsteps,(m*(j-1)+1):(m*j)] <- matrix(Y.update[(p+1-j):(nsteps+p-j),],ncol=m)
        }

      # Put on the exogenous regressors and make the constant the
      # first exog regressor after the AR coefs
      if(is.null(exog)==F)
        {
          X.update<-cbind(X.update,exog);
        }

      # Now, stack the original Y data and the augmented data.
      Y.update <- rbind(Y, ytmp[(capT+1):nrow(ytmp),])

      # Set up crossproducts and inverses we need
      X.update <- rbind(X, X.update)
      XX.update <- crossprod(X.update)     # Cross product of RHS variables
      hstar.update <- H0 + XX.update       # Prior + Cross product

      # Updated Regression estimates, Beta | Sigma
      B.update<-solve((hstar.update),(crossprod(X.update,Y.update) + H0[,1:m]))

      # Posterior mean of Sigma | Beta
      S.update <- (S0 + crossprod(Y.update)
                   + H0[1:m,1:m] - t(B.update)%*%(hstar.update)%*%(B.update))/(capT+nsteps+nu-m-1)

      # Posterior variance of B: VBh = diag(Sh.*.(inv((H0 + x'x)))
      hstarinv <- solve(hstar.update)
      vcv.Bh <- kronecker(S.update,hstarinv)

      # Draw from the conditional posterior pdfs of the parameters

      # This is only valid for just-identified models.
      df <- capT - m*p - m - 1 + nsteps
      wisharts <- rwishart(1, df, diag(m))

      # Generate the draws from the Wishart and the Beta
      # Wishart draw
      Sigmat <- (chol(S.update))
      Sigma.Draw <- t(Sigmat)%*%(df*solve(matrix(wisharts,m,m)))%*%Sigmat
      sqrtwish <- t(chol(Sigma.Draw))
      # Covariance of beta
      bcoefs.covar <- t(chol(vcv.Bh))

      # Draw of beta|Sigma ~ MVN(B.update, S.Update .*. Hstarinv)
      aplus <- matrix(B.update, ncol=1) +
        bcoefs.covar%*%matrix(rnorm(nrow(bcoefs.covar)), ncol=1)

      # Reshape and extract the coefs
      aplus <- matrix(aplus, ncol=m)
      aplus.coefs<-t(aplus[1:(m*p),]);          # extract the ar coefficients
      dim(aplus.coefs)<-c(m,m,p)                # push ar coefs into M x M x P array
      aplus.coefs<-aperm(aplus.coefs,c(2,1,3))  # reorder array

      intercept <- aplus[(m*p+1),]      # get drawn intercept....
      ar.coefs<-aplus.coefs             # AR coefs
      A0 <- sqrtwish

#      exog.coefs <- aplus[(m*p+2):nrow(aplus),]

      # Need to add something here to deal with the exogenous
      # regressors!


      # Print some intermediate results to capture progress....
      # and tell us that things are still running
      if (i%%500==0)
        { cat("Gibbs Iteration = ", i, "     \n");
          if(i<=burnin)
            { cat("(Still a burn-in draw.)\n");
            }


        }
      # Back to the top of the Gibbs loop....
    }
  endtime<-date()
  # Print time stamp so we know how long everything took.
  cat("Start time : ", starttime, "\n");
  cat("End time   : ", endtime, "\n");
    # Returns a list object
    output <- list(forecast=yforc, orig.y=y) #llf=ts(llf),hyperp=c(mu,prior)))
    attr(output, "class") <- c("forecast.VAR")
    return(output)
}

"plot.forecast.VAR" <-
function(x,y=NULL,varnames=NULL,
                               start=c(0,1),
                               freq=1, probs=c(0.05,0.95),
                               compare.level=NULL, ylab=NULL, ...)
{
    fcasts1 <- x
    fcasts2 <- y
    # compute quantities for ecdf of forecast matrix 1
    m <- dim(fcasts1$forecast)[3]
    h <- dim(fcasts1$forecast)[2]
    iters <- dim(fcasts1$forecast)[1]
    fcast1.summary <- array(apply(fcasts1$forecast, 3, forc.ecdf, probs=probs), c(h,3,m))

    # Now do the same for forecast 2 if non-NULL
    if (is.null(fcasts2)==FALSE)
      { fcast2.summary <- array(apply(fcasts2$forecast, 3,
                                      forc.ecdf, probs=probs),
                                c(h,3,m))
      }

#     # Now do the same for forecast 3 if non-NULL
#     if (is.null(fcasts3)==FALSE)
#       { fcast3.summary <- array(apply(fcasts3$forecast, 3,
#                                       forc.ecdf, probs=probs),
#                                 c(h,3,m))
#       }


    par(las=1, mar=c(1,2,2.5,1))
    for(i in 1:m)
    {
        forc1.ci <- ts(fcast1.summary[,,i], start=start)

        if(is.null(fcasts2)==FALSE){
            forc2.ci <- ts(fcast2.summary[,,i], start=start)
            forc.list <- c("forc1.ci","forc2.ci")
        } else {
            forc2.ci <- NULL
            forc.list <- c("forc1.ci")
        }

#         if(is.null(fcasts3)==FALSE)
#         { forc3.ci <- ts(fcast3.summary[,,i], start=start) }

        ylim <- c(floor(min(c(forc1.ci,forc2.ci,compare.level[i]))),
                  ceiling(max(c(forc1.ci,forc2.ci,compare.level[i]))))

        if(length(forc.list)==1){
            ts.plot(forc1.ci, gpars=list(lty=c(1,2,2), ylim=ylim, xlab="",axes=FALSE, ... ))
        } else if(length(forc.list==2)){
            ts.plot(forc1.ci, forc2.ci,
                    gpars=list(lty=c(1,1,1,2,2,2), ylim=ylim, xlab="",axes=FALSE, ... ))
        }

        axis(2,c(floor(min(c(forc1.ci,forc2.ci))), ceiling(max(c(forc1.ci,forc2.ci)))))
        mtext(varnames[i],side=3,line=1)

        box();
        if(i==1) { mtext(ylab, side=2, line=3, at=c(1.5*mean(ylim))) }
        abline(h=0)

        # put in the comparison level if one is provided
        if (is.null(compare.level)==FALSE)
            { abline(h=compare.level[i], lty=c(2)) }

      }
#    par(oldpar)
  }

