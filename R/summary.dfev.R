"summary.dfev" <-
function(object, latex=F, file=NULL, ...)
  { dfev.obj <- object
    errors <- dfev.obj$errors
    names <- dfev.obj$names
    std.err <- dfev.obj$std.err
    k <- dim(errors)[1]
    m <- dim(errors)[2]

    if(latex==T)
      { for (i in 1:m)
          { tmp <- matrix(errors[,,i], nrow=k, ncol=m)
            tmp <- cbind(std.err[,i],tmp)
            colnames(tmp) <- c("Std. Error", names)
            if(i==1)
              { # Ensure we clobber any old file
                print(xtable(tmp,
                             digits=rep(1,ncol(tmp)+1), 
                             caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                  file=file, append=F, table.placement="p")
              }
            else
              { print(xtable(tmp,
                         digits=rep(1,ncol(tmp)+1), 
                         caption=paste("Decomposition of Forecast Errors for a Shock to", names[i])),
                  file=file, append=T, table.placement="p")
              }
          }
      }
    else
      { for (i in 1:m)
          { cat(paste("Decomposition of Forecast Errors for a Shock to", names[i], "\n"))
                  cat("-------------------------------------------------------------\n")
            tmp <- matrix(errors[,,i], nrow=k, ncol=m)
            tmp <- cbind(std.err[,i],tmp)
            colnames(tmp) <- c("Std. Error", names)
            print(tmp)
                  cat("-------------------------------------------------------------\n")
          }
      }
  }

