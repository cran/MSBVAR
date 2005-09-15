"hc.forecast.var" <-
function(var.obj, yconst, nsteps,
                            burnin, gibbs, 
                            exog=NULL)
  {
    # Extract all the elements from the VAR object that we will need
    y <- var.obj$y
    ar.coefs <- var.obj$ar.coefs
    intercept <- var.obj$intercept
    exog.coefs <- var.obj$exog.coefs
    A0 <- t(chol(var.obj$mean.S))
    mu <- var.obj$hyperp
    prior <- var.obj$prior
    ncoef <- nrow(var.obj$Bhat)
    X <- var.obj$X         # rhs variables for the model
    Y <- var.obj$Y         # lhs variables for the model
    H0 <- var.obj$H0       # precision for the var coefs
    S0 <- var.obj$S0       # precision for the Sigma
    nu <- var.obj$prior[8] # df
    
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
      ytmp<-as.matrix(coef.forecast.var(y, intercept, ar.coefs, exog.coefs, m, p,
                                     capT, nsteps, A0=(A0)))
      
      # Get the impulse responses that correspond to the forecasted data
      M <- irf.var(var.obj, nsteps, A0)$mhat
      
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
  list(forecast=yforc, orig.y=y) #llf=ts(llf), hyperp=c(mu,prior)))
}

