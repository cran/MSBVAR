"forecast.var" <-
function(var.obj, nsteps, A0=t(chol(var.obj$mean.S)), 
                       shocks=matrix(0,nrow=nsteps,ncol=dim(var.obj$ar.coefs)[1]),
                       exog.fut=matrix(0,nrow=nsteps, ncol=nrow(exog.coefs)))
{
   # Set up the initial parameters for the VAR forecast function from
   #  VAR object
  y <- var.obj$y
  intercept <- var.obj$intercept
  ar.coefs <- var.obj$ar.coefs
  exog.coefs <- var.obj$exog.coefs
  m<-dim(ar.coefs)[1]
  p<-dim(ar.coefs)[3]
  capT<-nrow(y)
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

