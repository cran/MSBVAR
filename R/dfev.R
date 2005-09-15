"dfev" <-
function(var.obj, A0=t(chol(var.obj$mean.S)), k)
  { m <- dim(var.obj$ar.coefs)[1]
    p <- dim(var.obj$ar.coefs)[3]

    # Compute the IRF
    impulses <- irf.var(var.obj, k, A0=A0)$mhat
 
    # Find the cumulative innovations.  This is done by permuting the
    # array of the responses.  The irf function gives back an array
    # that m x m x k where k is the number of responses.  The last two
    # dimensions need to be permuted so we have the k x m x m arrays
    # of the responses.
    # The responses are also squared -- since they are in theory mean
    # zero, this gets them onto the right scale.
    impulses <- apply(aperm(impulses^2), c(2,3), cumsum)
    
    # Then compute the variances in each period.  This means we just
    # need to sum across the rows (forecast periods) for each shock
    # (now the outer array).
    
    var.imp <- apply(impulses,  c(1,3), sum)

    # Standardize the responses.
    for (i in 1:m)
      { impulses[,,i] <- impulses[,,i]/var.imp[,i] }
    
    # Scale into percentages
    errors <- 100*impulses
    # Output object
    output <- list(errors=errors, std.err=sqrt(var.imp),
                   names=colnames(var.obj$y))
    class(output) <- c("dfev")
    return(output)
  }

