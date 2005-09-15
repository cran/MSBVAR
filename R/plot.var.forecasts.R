"plot.var.forecasts" <-
function(x,y=NULL,varnames=NULL,
                               start=c(0,1), 
                               freq=1, probs=c(0.05,0.95),
                               compare.level=NULL, ylab=NULL, ...)
  {
    fcasts1 <- x
    fcasts2 <- y
    # compute quantities for ecdf of forecast matrix 1
    m <- dim(fcasts1$forecast)[3]
    h <- dim(fcasts1$forecast)[2]
    iters <- dim(fcasts1$forecast)[1]
    fcast1.summary <- array(apply(fcasts1$forecast, 3, forc.ecdf, probs=probs), c(h,3,m))

    # Now do the same for forecast 2 if non-NULL
    if (is.null(fcasts2)==FALSE)
      { fcast2.summary <- array(apply(fcasts2$forecast, 3,
                                      forc.ecdf, probs=probs),
                                c(h,3,m))
      }
    
#     # Now do the same for forecast 3 if non-NULL
#     if (is.null(fcasts3)==FALSE)
#       { fcast3.summary <- array(apply(fcasts3$forecast, 3,
#                                       forc.ecdf, probs=probs),
#                                 c(h,3,m))
#       }
    
    
    par(las=1, mar=c(1,2,2.5,1))
    for(i in 1:m)
      { forc1.ci <- ts(fcast1.summary[,,i], start=start)
        
        if(is.null(fcasts2)==FALSE)
          { forc2.ci <- ts(fcast2.summary[,,i], start=start) }

#         if(is.null(fcasts3)==FALSE)
#         { forc3.ci <- ts(fcast3.summary[,,i], start=start) }
        
        ylim <- c(floor(min(c(forc1.ci,forc2.ci,compare.level[i]))),
                  ceiling(max(c(forc1.ci,forc2.ci,compare.level[i]))))
        ts.plot(forc1.ci, forc2.ci, 
                 gpars=list(lty=c(1,1,1,2,2,2),
                            ylim=ylim, xlab="",axes=FALSE, ... ))
        axis(2,c(floor(min(c(forc1.ci,forc2.ci))),
                 ceiling(max(c(forc1.ci,forc2.ci)))))
        mtext(varnames[i],side=3,line=1)

        box();
        if(i==1) { mtext(ylab, side=2, line=3, at=c(1.5*mean(ylim))) }
        abline(h=0)

        # put in the comparison level if one is provided
        if (is.null(compare.level)==FALSE)
            { abline(h=compare.level[i], lty=c(2)) }
       
      }
#    par(oldpar)
  }

