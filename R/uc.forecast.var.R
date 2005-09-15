"uc.forecast.var" <-
function(var.obj, nsteps, burnin, gibbs,
                            exog=NULL)
  { # Extract all the elements from the VAR object
    y <- var.obj$y
    ar.coefs <- var.obj$ar.coefs
    intercept <- var.obj$intercept
    A0 <- t(chol(var.obj$mean.S))
    X <- var.obj$X         # rhs variables for the model
    Y <- var.obj$Y         # lhs variables for the model
    H0 <- var.obj$H0
    S0 <- var.obj$S0
    mu <- var.obj$hyperp
    exog.coefs <- var.obj$exog.coefs
    z <- var.obj$z
    lambda0 <- var.obj$prior[1]
    lambda1 <- var.obj$prior[2]
    lambda3 <- var.obj$prior[3]
    lambda4 <- var.obj$prior[4]
    lambda5 <- var.obj$prior[5]
    mu5 <- var.obj$prior[6]
    mu6 <- var.obj$prior[7]
    nu <- var.obj$prior[8]
    prior <- var.obj$prior.type
    num.exog <- var.obj$num.exog
    qm <- var.obj$qm
    ncoef <- nrow(var.obj$Bhat)
#     # build the extended regressions if necessary
#     if(is.na(sum(exog.fut))==F)
#       {
#         z.fut <- rbind(z,exog.fut)
#       }
#     else
#       { z.fut <- NULL
#       }
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
      ytmp <- forecast.var(var.obj, nsteps, A0=A0, shocks=epsilon.i)

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
  (list(forecast=yforc, orig.y=y, hyperp=c(mu,prior)))
}

