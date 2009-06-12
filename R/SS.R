## SS.R:  R-side functions for the SS C++ compressed state space class
## Date:  2007-11-15 -- Initial version
##        2008-10-14 -- Revised to use the bit package for storage.
##
## Description:
##
## Functions to draw samples of the state space for a model with filter
## probabilities 'fp'.  State spaces computed and stored in compressed
## form to reduce the computation and memory overhead involved in
## computing and storing state space draws.  R side storage of these
## state-spaces is done using the bit package.


# Draws a Markov state-space for a Markov-switching model.  This
# function does it all: filters, smoothes and samples the 0-1 elements
# of the state-space.
SS.draw <- function(xi, b, Q, m, h, n0, init.model, fp)
{
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")
    TT <- nrow(fp)
    SSe <- matrix(0, TT*m, h); Xik <- matrix(0, m*m, h)

    xi.tmp <- matrix(xi^2, m, h);
    for(i in 1:h){
        Xik[,i] <- as.vector(diag(1/xi.tmp[,i]))
        SSe[,i] <- as.vector(SS$e[,,i])
    }

    SQ <- steady.Q(Q)
    # Call C side SS draw function
    out <- .Call("SS.draw.cpp", SSe, Xik, Q, SQ, as.integer(TT),
                 as.integer(h), as.integer(m))
    return(out);
}

# Summary functions for MS state-spaces.

sum.SS <- function(x, ...)
{
    h <- x$h
    N <- length(x$ss.sample)
    ss <- x$ss.sample
    TT <- attributes(x$ss.sample[[1]])$n
    TTh <- TT*(h-1)

    # Now find the sums
    sums <- apply(matrix(unlist(lapply(x$ss.sample, as.integer)), nrow=TTh, ncol=N), 1, sum)
    sums <- matrix(sums, TT, h-1)
    sums <- cbind(sums, rep(N, TT) - rowSums(sums))  # Recover the h'th regime
    return(sums)
}

sum.SS1 <- function(x, c, ...)
{
    h <- x$h
    N <- length(x$ss.sample)
    ss <- x$ss.sample
    TT <- attributes(x$ss.sample[[1]])$n
    TTh <- TT*(h-1)
    sums <- matrix(0, TT, h)

    # Now find the sums
    tmp <- matrix(unlist(lapply(x$ss.sample, as.integer)), nrow=TTh,
                  ncol=N)


    for(i in 1:N)
    {
        tmp2 <- matrix(0, TT, h)
        tmp3 <- -1*(tmp[,i]-1)
        tmp2[,seq(c$cluster[i], by=ifelse(c$cluster[i]>1, -1, 1),
                  length.out=2)] <- cbind(tmp[,i], tmp3)
        sums <- sums + tmp2
    }
    return(sums)
}

mean.SS <- function(x, ...){
    sums <- sum.SS(x)
    N <- length(x$ss.sample)
    return(sums/N)
}


plot.SS <- function(x, ylab="State Probabilities", ...)
{
    tmp <- mean.SS(x)
    plot(ts(tmp), plot.type="single", col=1:ncol(tmp),
         ylim=c(0,1), ylab=ylab, ...)
    abline(h=0.5, lty=2, ...)

}


# Function to generate the vector of binary indicators of the state
# space.  This is a state-space generation function for one
# observation.
#
# p = filtered probabilities of each regime for an observation
# m = number of regimes.
#
# Returns a vector of 0-1 that indicate the regime.  A one is returned
# for the regime that the observation falls into.

## bingen <- function(p, Q, st1)
## {
##     h <- dim(Q)[1]
##     i <- 1

##     while(i<h)
##     { pr0 <- p[i]*Q[st1,i]/sum(p[i:h]*Q[st1,i:h])

##       if(runif(1)<=pr0)
##       { return(diag(h)[i,]) } else { st1 <- i <- i+1 }
##   }
##     return(diag(h)[h,])
##   }

## # Multi-move Gibbs sampler for the state space
## # filtered.prob  = BHLK.filter probabilities from BHLK.filter
## # P = Markov transition matrix

## generate.states <- function(filtered.prob, Q)
##   { TT <- nrow(filtered.prob)
##     h <- ncol(filtered.prob)

##     # storage
##     ss <- matrix(0, nrow=h, ncol=TT)

##     # generate the TT state
##     ss[,TT] <- bingen(filtered.prob[TT,], Q, 1)

##     for (t in (TT-1):1)
##       {
##           ss[,t] <- bingen(filtered.prob[t,], Q, which(ss[,t+1]==1))
##       }

##     return(t(ss))
##   }


## SS.draw.R <- function(xi, b, Q, m, h, n0, init.model, fp)
## {

##     SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

##     # Form Xi for the regimes
##     xi.tmp <- matrix(xi^2, m, h)
##     Xik <- array(0, c(m,m,h))

##     for(i in 1:h)
##     {

##         Xik[,,i] <- solve(diag(xi.tmp[,i]))
##     }

##     fp <- BHLK.filter(SS$e, Xik, Q)

##     ss <- generate.states(fp, Q)
##     return(ss)
## }

