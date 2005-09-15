"rwishart" <-
function(N, df, Sigma)
  { p = nrow(Sigma)
    SqrtSigma <- t(chol(Sigma))
    tmp <- array(0, c(p,p,N))
    for (i in 1:N)
      {
        Z <- matrix(rnorm(df*p),nrow=p,ncol=df)
        ZS <- crossprod(Z,SqrtSigma)
        tmp[,,i] <- crossprod(ZS)
      }
    if(N==1)
      {
        return(matrix(tmp,p,p))
      }
    else
      {
        return(tmp)
      } 
  }

