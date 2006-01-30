"irf.var" <-
function(var.obj, nsteps, A0=chol(var.obj$mean.S))
  {  
    ar.coef <- var.obj$ar.coef
    m<-dim(ar.coef)[1]                   # Capture the number of variables
     p<-dim(ar.coef)[3]                   # Capture the number of lags
     B<-array(0,c(m,m,max(nsteps,p)))     # Make an array for the AR coefs

     for(l in 1:(min(p,nsteps)))
       { B[,,l]<-ar.coef[,,l] }

     mhat<-matrix(0,ncol=m*m,nrow=nsteps) 
     dim(mhat)<-c(m,m,nsteps)             # Make an array to hold IRF
     mhat[,,1]<-A0                        # Identification condition
                                          # for IRF
     for(i in 2:nsteps)                   # Compute the IRF
       { for(j in 1:(i-1))
           { mhat[,,i]<-mhat[,,i] + (mhat[,,(i-j)]%*%B[,,j]) }
       }
    output <- list(B=B,mhat=mhat)
    class(output) <- c("irf.var")
    return(output)
  }

