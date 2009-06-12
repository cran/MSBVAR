
# Compute the steady state values (needed for the t=0 case) for a
# Markov process.  Not used here, since the initial state has to have
# probability 1/h

steady.Q <- function(P)
  { M <- dim(P)[1]
    if (M<3)
      {
        eta <- solve(rbind(cbind(1 - P[1,1], P[2,1]), rep(1,2)))%*%matrix(c(0,1))
      }
    else
      {
        eta <- solve(rbind(cbind(diag(M-1) - t(P)[1:(M-1),1:(M-1)], t(P)[1:(M-1),M]), rep(1,M)))%*%matrix(c(rep(0,M-1),1))
      }
    # Find the steady state if there are negative values -- that is
    # iterate a bit!
    while(cumprod(eta)[M]<0)
      {
        eta <- t(P)%*%eta
      }
    return(eta)
  }

####################################################################
# BHLK filter functions
# u = residuals for the regression model in each regime (T x m x h)
# Sigma = covariance of the residuals (m x m)
# Q = transition probability matrix (h x h)

BHLK.filter <- function(u, Sigma, Q)
{
    TT <- nrow(u); h <- nrow(Q); m <- dim(u)[2]
    u.vec <- matrix(0,TT*m,h); sig.vec <- matrix(0,m*m,h)
    for(i in 1:h){
        u.vec[,i] <- as.vector(t(u[,,i]))
        sig.vec[,i] <- as.vector(Sigma[,,i])
    }
    SQ <- steady.Q(Q)

    .Call("BHLK.filter.cpp", u.vec, sig.vec, Q, SQ, TT, h, m)
}

BHLK.filter.R <- function(u, Sigma, Q)
  { # Setup cconstants
    TT <- nrow(u)
    h <- nrow(Q)
    m <- dim(u)[2]
    pYSt <- matrix(NA, TT, h)
    filter.prob <- matrix(NA, nrow=h, ncol=TT)
    pSt <-  Q %*%steady.Q(Q)

    # Compute the likelihood for each observation in each regime
    for(i in 1:h) pYSt[,i] <- dmvnorm(u[,,i], sigma=Sigma[,,i])
    pYSt <- ifelse(pYSt==0, 1e-300, pYSt)

    for(t in 1:TT){
        num <- pYSt[t,]*pSt
        den <- pYSt[t,] %*% pSt
        filter.prob[,t] <- num/rep(den,h)
        # Update pSt
        pSt <- Q%*% filter.prob[,t]
    }

    return(t(filter.prob))
}
