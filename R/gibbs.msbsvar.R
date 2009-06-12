# gibbs.msbsvar.R
#
# Gibbs-MH sampler for MSBSVAR models
#
#
# Patrick T. Brandt
# 20081017 : Initial version

"gibbs.msbsvar" <- function(x, N1=1000, N2=1000, tune=matrix(10, x$m, x$h))
{

    # Initializations
    b <- x$b
    xi <- x$xi
    Q <- x$Q
    fp <- x$fp
    m <- x$m
    p <- x$p
    h <- x$h
    init.model <- x$init.model
    n0 <- x$n0
    n0cum <- x$n0cum
    abar <- x$abar
    bbar <- x$bbar
    alpha.prior <- x$alpha.prior

    # Modal A0(h) values around which to conduct likelihood
    # normalizations for the systems.

    mode <- array(0, c(m, m, h))

    b.tmp <- matrix(b, max(cumsum(n0)), h)

    for (i in 1:h)
    {
        mode[,,i] <- b2a(b.tmp[,i], init.model$Ui)
    }
    b.tmp <- matrix(NA, N1, length(b)+1)

# Burn-in

    Xik <- array(0, c(m,m,h));
    xi <- xi.tmp <- matrix(1, m, h)
    for(k in 1:h){ Xik[,,k] <- (diag(1/xi.tmp[,k]))}

    for (i in 1:N1)
    {
    # Update filter
        ss <- SS.draw(xi, b, Q, m, h, n0, init.model, fp)

        Q <- Q.draw(ss$transitions, prior=alpha.prior)

#    Xik <- array(0, c(m,m,h));  xi.tmp <- matrix(xi, m, h)
#    xi <- xi.draw(xi, b, m, h, n0, init.model, fp, abar, bbar)

        b <- b.draw(b, tune, m, h, n0, xi, init.model, fp, mode=mode, n0cum)$b
        b.tmp[i,] <- c(b, b.posterior(b, h, xi, init.model, fp, m, n0, n0cum))

        SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

        fp <- BHLK.filter(SS$e, Xik, Q)
        fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize

#        plot(ts(fp), plot.type="single", col=1:h, ylim=c(0,1))
#        abline(h=0.5)

        if(i%%1000==0) cat("Burnin iteration : ", i, "\n")
}

# normalized A0 at a peak
    b <- b.tmp[order(b.tmp[,ncol(b.tmp)],
                     decreasing=TRUE),][1,1:length(b)]

    A0init <- array(0, c(m, m, h))
    b.tmp <- matrix(b, max(cumsum(n0)), h)

    for (i in 1:h)
    {
        A0init[,,i] <- b2a(b.tmp[,i], init.model$Ui)
    }

# Storage matrices
    b.storage <- matrix(NA, N2, sum(n0)*h)
    F.storage <- matrix(NA, N2,  (h*m*(m*p + 1)))
    xi.storage <- matrix(1, N2, m*h)
    Q.storage <- matrix(NA, N2, h^2)
    accept.storage <- integer(length=m*h)
    ss.storage <- vector("list", N2)

    for (i in 1:N2)
    {
    # Update filter
        SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")
        fp <- BHLK.filter(SS$e, Xik, Q)
        fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]

    # Draw and save the state space
        ss <- SS.draw(xi, b, Q, m, h, n0, init.model, fp)

     # Compress and store the state-space
        ss.storage[[i]] <- as.bit.integer(ss$SS[,1:(h-1)])

    # Draw and save Q
        Q <- Q.draw(ss$transitions, prior=alpha.prior)
        Q.storage[i,] <- as.vector(Q)

    # Draw and save b
        b.tmp <- b.draw(b, tune, m, h, n0, xi, init.model, fp,
                        mode=A0init, n0cum)
        b.storage[i,] <- b <- b.tmp$b
        accept.storage <- as.integer(accept.storage + b.tmp$accept)

    # Compute and save F
        SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")
        F.storage[i,] <- as.vector(SS$F)

    # Draw and save xi
#    xi.storage[i,] <- xi <- xi.draw(xi, b, m, h, n0, init.model, fp, abar, bbar)
#    Xik <- array(0, c(m,m,h));  xi.tmp <- matrix(xi, m, h)
#    for(k in 1:h){ Xik[,,k] <- (diag(1/xi.tmp[,k]))}


        if(i%%1000==0) {
            cat("Final iteration : ", i, "\n")
#            plot(ts(fp), plot.type="single", col=1:h, ylim=c(0,1))
#            abline(h=0.5)

        }
    }


    # Now class things to we can later generate summary methods
    class(ss.storage) <- c("SS")

    output <- list(b.sample=mcmc(b.storage),
                   F.sample=mcmc(F.storage),
                   xi.sample=mcmc(xi.storage),
                   Q.sample=mcmc(Q.storage),
                   ss.sample=ss.storage,
                   accept.rate=accept.storage/N2,
                   A0mode=A0init,
                   h=h)

    class(output) <- c("MSBSVAR")

    return(output)

}


# Single draws of the model parameters

# b.draw -- draws one vector of b using a Metropolis-Hastings
# algorithm.  Saves the computations from the "selected" draw to speed
# computation in the next draw.

b.draw <- function(b, tune, m, h, n0, xi, init.model, fp, mode, n0cum)
{

    # get moments for "old"
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.b")

    b <- matrix(b, ncol=h)
    bout <- b
    acceptout <- matrix(NA, m, h)

#    df <- (colSums(round(fp)))+1

    A0k <- array(0, c(m,m,h))
    for(k in 1:h)
    { for (j in 1:m)
      {

          # Draws
          bcov <- tune[j,k]*(SS$bjk.cov[[k]][[j]])
          bold <- b[(n0cum[j]+1):(n0cum[(j+1)]),k]
          bnew <- rmvnorm(1, sigma=bcov)

          # Metropolis step
          p.old <- dmvnorm(bold, sigma=bcov)
          p.new <- dmvnorm(bnew, sigma=bcov)
          alpha <- min(1, p.new/p.old)

          if(runif(1)<alpha) {
              bout[(n0cum[j]+1):(n0cum[(j+1)]),k] <- bnew
              acceptout[j,k] <- 1
          } else {  acceptout[j,k] <- 0  }
      }
      # now normalize the draw
      A0k[,,k] <- b2a(bout[,k], init.model$Ui)
      A0k[,,k] <- normalize.svar(A0k[,,k], mode[,,k], method="DistanceMLA")$A0normalized
      bout[,k] <- a2b(A0k[,,k], init.model$Ui)
  }


    return(list(b=as.vector(bout), accept=as.vector(acceptout)))
}


## b.draw.init <- function(b, tune, m, h, n0, xi, init.model, fp)
## {

##     # get moments for "old"
##     SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.b")

##     n0cum <- c(0, cumsum(n0))
##     b <- matrix(b, ncol=h)
##     bout <- b
##     acceptout <- matrix(NA, m, h)
##     df <- (colSums(round(fp)))+1
##     A0k <- array(0, c(m,m,h))

##     for(k in 1:h)
##     { for (j in 1:m)
##       {

##           # Draws
##           bcov <- tune[j,k]*solve(SS$bjk.cov[[k]][[j]])/df[k]
##           bold <- b[(n0cum[j]+1):(n0cum[(j+1)]),k]
##           bnew <- rmvnorm(1, sigma=bcov)

##           # Metropolis step
##           p.old <- dmvnorm(bold, sigma=bcov)
##           p.new <- dmvnorm(bnew, sigma=bcov)
##           alpha <- min(1, p.new/p.old)

##           if(runif(1)<alpha) {
##               bout[(n0cum[j]+1):(n0cum[(j+1)]),k] <- bnew
##               acceptout[j,k] <- 1
##           } else {  acceptout[j,k] <- 0  }
##       }

##   }

##     return(list(b=as.vector(bout), accept=as.vector(acceptout)))
## }


xi.draw <- function(xi, b, m, h, n0, init.model, fp, abar, bbar)
{
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.xi")

    xi <- matrix(0, m, h)

    for (i in 1:h)
    {
        for (j in 1:m)
        {
            fjk <- SS$F[,j,i]
            alpha <- abar[j,i] + 0.5*SS$moments[[i]][[1]]
##             beta <- bbar[j,i] + 0.5*((crossprod(SS$A0k[,j,i], SS$moments[[i]][[2]])/xi[j,i] %*%
##                                       SS$A0k[,j,i]) -
##                                      2*crossprod(fjk,SS$moments[[i]][[3]])%*%matrix(SS$A0k[,j,i], ncol=1)/xi[j,i]
##                                      + crossprod(fjk,
##                                                  SS$moments[[i]][[4]])%*%fjk)/xi[j,i]


            ssr <- sum((SS$moments[[i]][[5]]%*%SS$A0k[,j,i] - SS$moments[[i]][[6]]%*%SS$F[,j,i])^2)
            beta <- bbar[j,i] + 0.5*ssr

            xi[j, i] <- rgamma(1, shape=alpha, scale=1/beta)
        }

    }
    # normalize the xi
    xi <- norm.xi(xi)

    return(as.vector(xi))

}



Q.draw <- function(SStrans, prior)
{
    # compute the density based on the counts.
    h <- ncol(SStrans)
    alpha <- SStrans + prior - matrix(1, h, h)
    Q <- t(apply(alpha, 1, rdirichlet, n=1))
    return(Q)
}



