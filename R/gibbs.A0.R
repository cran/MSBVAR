"gibbs.A0" <- function(varobj, N1, N2, thin=1, normalization="DistanceMLA")
{
    # Leave the character strings on the R side and pass an int flag
    # to the normalization routine norm_svar (fully documented in
    # MSBVARfun.cpp) in the order below

    tmp <- sanity.check.gibbs(list(N1=N1, N2=N2, thin=thin, normalization=normalization))
    methodlist <- c("DistanceMLA", "DistanceMLAhat", "Euclidean", "PositiveDiagA", "PositiveDiagAinv")

    if(tmp){ method <- which(methodlist==normalization)-1 } else {method <- which(methodList==tmp)-1}

    cat("Normalization Method: ", normalization, "(", method,")\n")

    .Call("gibbsA0.cpp", varobj, as.integer(N1), as.integer(N2),
          as.integer(thin), as.integer(method), gibbs.setup.bsvar(varobj)$UT)

}


"A02mcmc" <- function(x)
{
    return(mcmc(matrix(x$A0.posterior$A0,
                       nrow=x$N2,
                       ncol=length(x$A0.posterior$struct),
                       byrow=T),
                thin=x$thin))
}


## "gibbs.A0.DEPRECATED" <- function(varobj, N1, N2, thin=1, normalization="DistanceMLA")
## {return(gibbs.A0.BSVAR(varobj, N1, N2, thin=1, normalization="DistanceMLA"))}

## # Gibbs sampler for A0's in B-SVAR model

## "gibbs.A0.BSVAR.DEPRECATED" <- function(varobj, N1, N2, thin=1, normalization="DistanceMLA")
## {
##     A0ml <- varobj$A0.mode
##     setup <- gibbs.setup.bsvar(varobj)
##     UT <- setup$UT
##     df <- varobj$df
##     n0 <- varobj$n0
##     ident <- varobj$ident
##     m <- ncol(A0ml)
##     a0indx <- which(ident==0)  # index of the zeros in A0

##     # storage for the A0.posterior
##     A0.flat <- vector(mode="numeric",length=(sum(n0)*N2))
##     A0.struct <- which(ident!=0)
##     A0.posterior <- list(A0=A0.flat,struct=A0.struct,m=m)

##     W.flat <- vector(mode="numeric", length=((m-1)*(m-1)*m*N2))
##     W.index <- vector(mode="numeric", length=(m*N2))

##     # counter for thinned draws
##     keep <- 1

##     # MAIN GIBBS SAMPLER LOOP FOR SVAR MODELS

##     # Normalize the A0 ML estimate to have a positive diagonal --
##     # easier for interpretation

##     A0ml <- normalize.svar(varobj$A0.mode, varobj$A0.mode,
##                            method=normalization)$A0normalized
##     A0gbs <- A0ml

##     for (i in 1:N1)
##     {
##         W.tmp <- vector(mode="list",length=m)
##         drawi <- drawA0(A0gbs, UT, df, n0, W.tmp)
##         A0gbs <- drawi$A0gbs
##         W.tmp <- drawi$W

##         # clean up A0gbs so we minimize the roundoff error.
##         A0gbs[a0indx] <- 0

##         # Now normalize the A0 draw
##         A0normalized <- normalize.svar(A0gbs, A0ml,
##                                        method=normalization)

##         A0gbs <- A0normalized$A0normalized

##         if (i%%100==0)
##         { cat("Gibbs Burnin Iteration = ", i, "\n"); }
##     }

##     switch.counter <- 0

##     for (i in 1:(N2*thin))
##     {
##         drawi <- drawA0(A0gbs, UT, df, n0, W.tmp)
##         A0gbs <- drawi$A0gbs
##         W <- W.flatten(drawi$W)
##         W.tmp <- drawi$W

##       # clean up A0gbs so we minimize the roundoff error.
##         A0gbs[a0indx] <- 0

##       # Now normalize the A0 draw
##         A0normalized <- normalize.svar(A0gbs, A0ml,
##                                        method=normalization,
##                                        switch.count=switch.counter)
##       # Save the thin'th normalized A0
##         if(i%%thin==0)
##         {
##             if(keep==1)
##             {
##                 W.index[1:m] <- W$W.index
##                 W.flat[1:W.index[m]] <- W$W
##                 len <- length(A0.flatten(A0normalized$A0normalized))
##                 A0.posterior$A0[1:len] <- A0.flatten(A0normalized$A0normalized)
##             }
##             else
##             {
##                 st <- (keep-1)*m+1
##                 ed <- keep*m

##                 W.index[st:ed] <- W.index[st-1] + W$W.index
##                 W.flat[(W.index[st-1]+1):W.index[ed]] <- W$W

##                 st <- (keep-1)*length(A0.posterior$struct)
##                 ed <- st+length(A0.posterior$struct)
##                 A0.posterior$A0[(st+1):ed] <- A0.flatten(A0normalized$A0normalized)
##             }
##             keep <- keep + 1
##         }

##       # Reinitalize
##         A0gbs <- A0normalized$A0normalized
##       # Iterate switch counter
##         switch.counter <- switch.counter + A0normalized$switch.count

##         if (i%%100==0)
##         {
##           cat("Gibbs Iteration = ", i, "\n")
##           ld <- determinant(A0normalized$A0normalized)
##           cat("A0 log-det.     = ", ld$sign*ld$modulus, "\n")
##         }

##     }
##     W.posterior <- list(W=W.flat[1:(W.index[(keep-1)*m])], W.index=W.index, m=m)
##     return(list(A0.posterior=A0.posterior, W.posterior=W.posterior,
##                 ident=ident, thin=thin, N2=N2))
## }



# Old version without conserving storage of A0's and W's
## "gibbs.A0.structural.szbvar" <- function(varobj,
##                                          N1, N2, thin=1, normalization="PositiveDiagA")
##   {
##       A0ml <- varobj$A0.mode
##       setup <- gibbs.setup.bsvar(varobj)
##       UT <- setup$UT
##       df <- varobj$df
##       n0 <- varobj$n0
##       ident <- varobj$ident
##       m <- ncol(A0ml)
##       a0indx <- which(ident==0)  # index of the zeros in A0

##       # storage for the A0.posterior
##       A0.posterior <- array(0, c(m,m,N2))
##       W.posterior <- list("vector", N2)

##       # counter for thinned draws
##       keep <- 1

##       # MAIN GIBBS SAMPLER LOOP FOR SVAR MODELS

##       # Normalize the A0 ML estimate to have a positive diagonal --
##       # easier for interpretation

##       A0ml <- normalize.svar(varobj$A0.mode, varobj$A0.mode,
##                              method=normalization)$A0normalized
##       A0gbs <- A0ml

##       for (i in 1:N1)
##       {     W <- vector("list", m)
##             drawi <- drawA0(A0gbs, UT, df, n0, W)
##             A0gbs <- drawi$A0gbs
##             W <- drawi$W
##          # clean up A0gbs so we minimize the roundoff error.
##          A0gbs[a0indx] <- 0

##          # Now normalize the A0 draw
##          A0normalized <- normalize.svar(A0gbs, A0ml,
##                                         method=normalization)

##          A0gbs <- A0normalized$A0normalized

##          if (i%%100==0)
##            { cat("Gibbs Burnin Iteration = ", i, "\n"); }
##        }

##     switch.counter <- 0

##     for (i in 1:(N2*thin))
##       {     drawi <- drawA0(A0gbs, UT, df, n0, W)
##             A0gbs <- drawi$A0gbs
##             W <- drawi$W

##       # clean up A0gbs so we minimize the roundoff error.
##             A0gbs[a0indx] <- 0

##       # Now normalize the A0 draw
##             A0normalized <- normalize.svar(A0gbs, A0ml,
##                                            method=normalization,
##                                            switch.count=switch.counter)
##       # Save the thin'th normalized A0
##             if(i%%thin==0)
##               {
##                 A0.posterior[,,keep] <- A0normalized$A0normalized
##                 W.posterior[[keep]] <- W
##                 keep <- keep + 1
##               }

##       # Reinitalize
##             A0gbs <- A0normalized$A0normalized
##       # Iterate switch counter
##             switch.counter <- switch.counter + A0normalized$switch.count

##             if (i%%100==0)
##               { cat("Gibbs Iteration = ", i, "\n")
##                 ld <- determinant(A0normalized$A0normalized)
##                 cat("A0 log-det.     = ", ld$sign*ld$modulus, "\n")
##               }

##           }
##     return(list(A0.posterior=A0.posterior, W.posterior=W.posterior,
##                 ident=ident, thin=thin, N2=N2))
##   }
