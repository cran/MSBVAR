"mc.irf.var" <-
function(varobj, nsteps,draws)
{ m<-dim(varobj$ar.coefs)[1]  # Capture the number of variablesbcoefs <- varobj$Bhat
  p<-dim(varobj$ar.coefs)[3]    # Capture the number of lags
  
  bcoefs <- t(varobj$Bhat[1:(m*p),])
  dim(bcoefs) <- c((m^2*p),1)
  capT <- nrow(varobj$Y)
  X <- varobj$X[,1:(m*p)]
  XXinv <- solve(crossprod(X))

  # Get the correct moment matrix for the BVAR object
  if(is.null(varobj$prior)==FALSE)
    { XXinv <- solve(varobj$hstar[1:(m*p),1:(m*p)]) }

  # Cholesky of the covariance
  Sigmat <- chol(varobj$mean.S)
  
  # Matrices to hold stuff
  impulse <- matrix(0,nrow=draws,ncol=(m^2*nsteps))

  # Do all the Wishart draws
  # DF from Box and Tiao Section 8.5.1 p. 460. and Sims and Zha 1998.
  df <- capT - m*p - m - 1
  wisharts <- rwishart(draws, df, diag(m))
  XXinv <- t(chol(XXinv))
  
  # Main loop -- uses antithetic acceleration to cut the number of
  # effective draws in half.
  
  for(i in 1:draws)
    {
      # Generate the draws from the Wishart and the Beta
      Sigma.Draw <- t(Sigmat)%*%(capT*solve(matrix(wisharts[,,i],m,m)))%*%Sigmat
      sqrtwish <- t(chol(Sigma.Draw))
      
           # Here we exploit the fact that since there is a Kronecker
           # product structure to the covariance mtx of beta.  This means
           # we can take the Cholesky of each component in the Kronecker
           # product of the VCV matrix of Beta and combine them with random
           # normal draws to get a draw from beta.  This is much much easier
           # than using a draw from a full MVN pdf

          bcoefs.covar <- kronecker(sqrtwish, XXinv)
          coef.u <- bcoefs.covar%*%matrix(rnorm(m^2*p), ncol=1)
          coef.mtx.odd <- bcoefs + coef.u
          coef.mtx.even <- bcoefs - coef.u

          # Now compute the irfs
          impulse[i,] <- t(irf.var.from.beta(sqrtwish,
                                             coef.mtx.odd,
                                             nsteps))
      if (i%%100==0)
        { cat("Monte Carlo IRF Integration Iteration = ", i, "\n"); }
    }
  
  # Put the results into an array of the impulses.
  output <- array(impulse,c(draws,nsteps,m^2))
  class(output) <- c("mc.irf.var")
  return(output)
}

