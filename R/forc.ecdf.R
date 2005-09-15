"forc.ecdf" <-
function(forecasts, probs=c(0.05,0.95), start=c(0,1), ...)
  { m.forecast<-apply(forecasts,2,mean)
    quant<-apply(forecasts, 2, quantile, probs)
    vplus<-quant[1,]
    vminus<-quant[2,]
    forc.ci<-ts(t(rbind(m.forecast,vplus,vminus)),start=start, ...)
    return(forc.ci)
  }

