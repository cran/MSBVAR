"plot.irf.var" <-
function(x,varnames=NULL, ...)
  { impulses<-x$mhat               # Push IRF results into another array
    m<-dim(impulses)[1]              # Get array dimensions
    nsteps<-dim(impulses)[3]
    dim(impulses)<-c((m^2),nsteps)   # Reshape array to get a matrix of irf series

    # Now make this into time series
    impulses<-ts(t(impulses),start=c(0,1),freq=1)
    
    # Get the min and max for the rows of the irfs
    dim(x$mhat) <- c(m,(m*nsteps))
    minmax <- apply(x$mhat,1,range)

    # Graph the irf 
    par(mfrow=c(m,m),mai=c(0.25,0.25,0.15,0.25), omi=c(0.15,0.75,1,0.15))
    for (i in 1:m)
      { for (j in 0:(m-1))
          { col <- m*j + i
              plot(impulses[,col],
                 xlab="",ylab="",
                 ylim=c(minmax[,i])
                 )
            abline(h=0)
          # Add plot axis labels
          if(j==0)
            { mtext(varnames[i],side=2,line=3) }
          if(i==1)
            { mtext(varnames[j+1],side=3,line=2) }

      }
      }
    mtext("Shock to", side=3, line=3, outer=T)
    mtext("Response in", side=2, line=3,, outer=T)
#    return(impulses)
  }

