# Updated 2006-01-29 to transpose rows and columns so the responses
# are in the rows.

"plot.irf.var" <-
function(x,varnames=NULL, ...)
  { impulses<-x$mhat               # Push IRF results into another array
    m<-dim(impulses)[1]              # Get array dimensions
    nsteps<-dim(impulses)[3]
    stacked.impulses <- impulses
    dim(stacked.impulses) <- c(m*nsteps,m)
    dim(impulses)<-c((m^2),nsteps)   # Reshape array to get a matrix of irf series

    # Now make this into time series
    impulses<-ts(t(impulses),start=c(0,1),freq=1)

    # Get the min and max for the rows of the irfs
    minmax <- apply(stacked.impulses, 2, range)

    # Counter for varnames
    j <- 1
    
    # Graph the irf
    par(mfrow=c(m,m),mai=c(0.25,0.25,0.15,0.25), omi=c(0.15,0.75,1,0.15))
    for (i in 1:m^2)
      { 
        plot(impulses[,i],
             xlab="",ylab="",
             ylim=c(minmax[,ceiling(i/m)])
             )
        abline(h=0)
        
        # Add plot axis labels
        if(i<=m)
          { mtext(varnames[i],side=3,line=2) }
        if((i-1)%%m==0)
          { mtext(varnames[j],side=2,line=3)
            j <- j+1
          }

      }
    mtext("Shock to", side=3, line=3, outer=T)
    mtext("Response in", side=2, line=3,, outer=T)
  }

