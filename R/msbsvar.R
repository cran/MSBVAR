# msbsvar.R
# MS-BSVAR model code and summary functions
#
# Patrick T. Brandt
#
# 20081017 : Initial version

# Mode finder for MSBSVAR

# Need a sanity.check() for this function!
"msbsvar" <- function(Y, z=NULL, p, h, ident,
                      lambda0, lambda1, lambda3,
                      lambda4, lambda5, mu5, mu6, qm,
                      alpha.prior,
                      max.iter=10)
{

    m <- ncol(Y)
    # Set up and fit initializing szbsvar() model
    init.model <- szbsvar(Y, z=z, p=p, lambda0=lambda0,
                          lambda1=lambda1, lambda3=lambda3,
                          lambda4=lambda4, lambda5=lambda5,
                          mu5=mu5, mu6=mu6, qm=qm,
                          ident=ident)
    # Get the residuals
    u <- init.model$structural.innovation[(m+1):nrow(init.model$structural.innovation),]
    u <- array(u, c(nrow(u), ncol(u), h))

    # Set up h-regime values of the model: A0(h), xi(h), Q, Sigma(h), b(h)

    Q <- (1-(h*0.01/(h-1)))*diag(h) + matrix(0.01/(h-1), h, h)
    A0 <- init.model$A0.mode
    xi <- rep(1, m*h)

    A0h <- array(0, c(m,m,h))
    for(i in 1:h) for(j in 1:m)  A0h[,j,i] <- A0[,j]

    Sigma <- array(apply(A0h, 3, crossprod), c(m,m,h))

    n0 <- init.model$n0

    abar <- bbar <- matrix(1, m, h)
    n0cum <- c(0, cumsum(init.model$n0))

    # Stack up the free parameters : will need to follow this convention
    # in all of the objective function inputs
    b <- rep(init.model$b, h) + rnorm(sum(ident)*h)

    # Inital filter

    fp <- rdirichlet(nrow(u), rep(1, h))

    # Find the peak of the posterior / mode.
    posterior.mode <- MSBSVAR.posterior.block(Y, h, p, max.iter=max.iter,
                                              init.model, b, xi, fp, Q, m,
                                              n0, n0cum, alpha.prior,
                                              abar, bbar)

    return(posterior.mode)

}





################################################################
# Other helper function / interim computations
################################################################

# Function to normalize xi vectors
norm.xi <- function(xi)
{
    return(xi/matrix(xi[1,], nrow(xi), ncol(xi), byrow=TRUE))

}

# Posterior functions (equation 18 and p. 23 of SWZ '06).  Break this
# into three optimands as described in SWZ, p. 27-28.

"b.posterior" <- function(b, h, xi, init.model, fp, m, n0, n0cum)
{
    fn <- 0.0
    b <- matrix(b, max(n0cum), h)

    # Get regime specific SS
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.xi")

    # Prior covars

    covar <- init.model$H0inv.tilde #, solve)

    for(i in 1:h)
    {
        for(j in 1:m)
        {

            fjk <- SS$F[,j,i]
            btmp <- b[(n0cum[j]+1):(n0cum[(j+1)]),i]

            # Get the log det we need.
            ld <- determinant(SS$A0k[,,i])
            ld <- ld$sign*ld$modulus

            fn <- fn + -0.5* crossprod(btmp, covar[[j]])%*%btmp + ld*SS$moments[[i]][[1]] +
                -0.5*(crossprod(SS$A0k[,j,i], SS$moments[[i]][[2]]) %*%
                      SS$A0k[,j,i]  -
                          2*crossprod(fjk,SS$moments[[i]][[3]])%*%matrix(SS$A0k[,j,i], ncol=1) +
                              crossprod(fjk, SS$moments[[i]][[4]])%*%fjk)

        }

    }


    return(fn)
}



# Note that the xi^2 are from the fractional version of the Gamma pdf,
# so the moments are E[xi^2] = a/b

xi.posterior <- function(xi, h, m, b, n0, init.model, fp, abar, bbar)
{

    # Get regime specific SS
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model,
                     method="update.xi")

    # normalize the xi
    xi <- norm.xi(matrix(xi, m, h))

    fn <- 0.0

    for (i in 1:h)
    {
        for (j in 1:m)
        {
            fjk <- SS$F[,j,i]
            alpha <- abar[j,i] + 0.5*SS$moments[[i]][[1]]
##             beta <- bbar[j,i] + 0.5*((crossprod(SS$A0k[,j,i], SS$moments[[i]][[2]]) %*%
##                                       SS$A0k[,j,i]) -
##                                      2*crossprod(fjk,SS$moments[[i]][[3]])%*%matrix(SS$A0k[,j,i], ncol=1)
##                                      + crossprod(fjk,
##                                                  SS$moments[[i]][[4]])%*%fjk)

            ssr <- sum((SS$moments[[i]][[5]]%*%SS$A0k[,j,i] - SS$moments[[i]][[6]]%*%SS$F[,j,i])^2)
            beta <- bbar[j,i] + 0.5*ssr
            fn <- fn + dgamma(xi[j,i], shape=alpha, scale=1/beta, log=TRUE)
        }
    }

    return(fn)
}

##########################################################
# Q.posterior : objective function or Pr(Q | S_t, Y)
# Optimization of this is done below.
##########################################################

Q.posterior <- function(q.init, xi, h, m=m, fp, b, n0, init.model, prior)
{
    n0cum <- c(0,cumsum(n0))
    b <- matrix(b, max(n0cum), h)

    # Get residuals
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

    # Get the structural residuals

    fn <- 0.0

    # normalize xi
    xi <- norm.xi(matrix(xi, m, h))

    Xik <- array(0, c(m,m,h))
    # Form A(k) for the regimes
    # Form Xi for the regimes

    for(i in 1:h)
    {
        Xik[,,i] <- (diag(1/xi[,i]))
    }

    # filter the data based on the earlier results
    fp <- BHLK.filter(SS$e, Xik, matrix(q.init, h, h))

    # count the transitions for the newly filtered data
    ct <- count.transitions(round(fp))

    # compute the density based on the counts.
    alpha <- ct + prior - matrix(1, h, h)
    fn <- sum(log(ddirichlet(matrix(q.init, h, h), alpha)))

    return(fn)

}

##########################################################
# Computes the gradient for Q.posterior()
##########################################################

Q.grad <- function(q.init, xi, h, m=m, fp, b, n0, init.model, prior)
{
    n0cum <- c(0,cumsum(n0))
    b <- matrix(b, max(n0cum), h)

    # Get residuals
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

    # Get the structural residuals

    fn <- 0.0
    xi <- norm.xi(matrix(xi, m, h))

    Xik <- array(0, c(m,m,h))
    # Form A(k) for the regimes
    # Form Xi for the regimes

    for(i in 1:h)
    {
        Xik[,,i] <- (diag(1/xi[,i]))
    }

    # filter the data based on the earlier results
    fp <- BHLK.filter(SS$e, Xik, matrix(q.init, h, h))

    # count the transitions for the newly filtered data
    ct <- count.transitions(round(fp))

    # compute the density based on the counts.
    alpha <- ct + prior - matrix(1, h, h)

    return(q.init/as.vector(alpha))

}

# Counts the transitions in the discrete state space.

count.transitions <- function(s)
  { M <- ncol(s)
    TT <- nrow(s)
    s <- crossprod(t(s), as.matrix(seq(1:M)))
    sw <- matrix(0, M, M)
    for (t in 2:TT)
      { st1 <- s[t-1]
        st <- s[t]
        sw[st1,st] <- sw[st1, st] + 1
      }
    return(sw)
  }



Q1.posterior <- function(q.init, xi, h=h, m, fp, b, n0, init.model, prior)
{
    Q.full <- cbind(matrix((q.init), h, h-1), matrix(1, h, 1) - rowSums(matrix(q.init, h, h-1)))
    n0cum <- c(0,cumsum(n0))
    b <- matrix(b, max(n0cum), h)

    # Get residuals
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

    # Get the structural residuals

    fn <- 0.0
    xi <- norm.xi(matrix(xi, m, h))

    Xik <- array(0, c(m,m,h))
    # Form A(k) for the regimes
    # Form Xi for the regimes

    for(i in 1:h)
    {
        Xik[,,i] <- (1/diag(xi[,i]))
    }

    # filter the data based on the earlier results
    fp <- BHLK.filter(SS$e, Xik, matrix(Q.full, h, h))
    fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize

    # count the transitions for the newly filtered data
    ct <- count.transitions(round(fp))

    # compute the density based on the counts.
    alpha <- ct + prior - matrix(1, h, h)
    fn <- sum(log(ddirichlet(matrix(Q.full, h, h), alpha)))

    return(fn)

}



Q1.grad <- function(q.init, xi, h=h, m=m, fp, b, n0, init.model, prior)
{
    Q.full <- cbind(matrix((q.init), h, h-1), matrix(1, h, 1) - rowSums(matrix(q.init, h, h-1)))
    n0cum <- c(0,cumsum(n0))
    b <- matrix(b, max(n0cum), h)

    # Get residuals
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

    # Get the structural residuals

    fn <- 0.0
    xi <- norm.xi(matrix(xi, m, h))

    Xik <- array(0, c(m,m,h))
    # Form A(k) for the regimes
    # Form Xi for the regimes

    for(i in 1:h)
    {
        Xik[,,i] <- (diag(1/xi[,i]))
    }

    # filter the data based on the earlier results
    fp <- BHLK.filter(SS$e, Xik, Q.full)

    # count the transitions for the newly filtered data
    ct <- count.transitions(round(fp))

    # compute the density based on the counts.
    alpha <- ct + prior - matrix(1, h, h)

    return(q.init/as.vector(alpha[,1:(h-1)]))

}


# Optimize the MS-B-SVAR model by blocks.

MSBSVAR.posterior.block <- function(Y, h, p, max.iter=10, init.model, b, xi,
                                    fp, Q,
                                    m, n0, n0cum,
                                    alpha.prior, abar, bbar)
{

    # Control parameters for the optimizations
    op.method <- c("Nelder-Mead", "BFGS", "L-BFGS-B")
    iter <- 1

    # Initial function values

    max.b.fn <- b.posterior(b, h, xi, init.model, fp, m, n0, n0cum)
#    max.xi.fn <- xi.posterior(xi, h, m, b, n0, init.model=model, fp, abar,
#                              bbar)
    max.Q.fn <- Q.posterior(as.vector(Q), xi, h, m=m, fp, b, n0,
                            init.model, prior=alpha.prior)

    n0cum <- c(0,cumsum(n0))

    for(iter in 1:max.iter)
    {

    b.old <- b
    xi.old <- xi
    Q.old <- Q

    # Optimize the Q block


     q.init <- as.vector(Q[,1:(h-1)])

    max.Q <- optim(par = q.init,
                   f = Q1.posterior,
                   gr = Q1.grad,
                   xi=xi, h=h,
                   m=m, fp=fp, b=b, n0=n0, init.model=init.model,
                   prior=alpha.prior,
                   method = "L-BFGS-B",
                   lower = rep(0.001, h),
                   upper = rep(0.999, h),
                   control=list(trace=1, fnscale=-nrow(Y),
                   lmm=10,
                   maxit=500))

    Q <- cbind(matrix(max.Q$par, h, h-1), matrix(1, h, 1) -
               rowSums(matrix(max.Q$par, h, h-1)))

    xi.tmp <- matrix(xi, m, h)

    Xik <- array(0, c(m,m,h))
    btmp <- matrix(b, max(n0cum), h)

    for(i in 1:h)
    {

        Xik[,,i] <- (diag(1/xi.tmp[,i]))
    }

     SS <- hregime.SS(h, m, fp, b, xi, n0, init.model=init.model,
                      method="update.Q")

     fp <- BHLK.filter(SS$e, Xik, Q)

     fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize



##     try(max.xi <- optim(xi, xi.posterior,
##                     method=op.method[2], hessian=FALSE,
##                     control=list(maxit=200, fnscale=-nrow(Y),
##                     trace=1, REPORT=20), #ifelse(iter<3, 0, 1)),
##                     h=h, m=m, b=b, n0=n0, init.model=model, fp=fp,
##                     abar=abar, bbar=bbar))

##     xi <- norm.xi(matrix(max.xi$par, m, h))


    try(max.b <- optim(b, b.posterior,
                   method=op.method[2], hessian=FALSE,
                   control=list(maxit=200, fnscale=-nrow(Y),
                   trace=1, REPORT=20), #ifelse(iter<3, 0, 1)),
                   h=h, xi=xi, init.model=init.model,
                   fp=fp, m=m, n0=n0, n0cum=n0cum))

    b <- max.b$par

##     # Form Xi for the regimes
##     xi.tmp <- xi

##     Xik <- array(0, c(m,m,h))
##     btmp <- matrix(b, max(n0cum), h)

##     for(i in 1:h)
##     {

##         Xik[,,i] <- (diag(1/xi.tmp[,i]))
##     }

    # Update filter
    SS <- hregime.SS(h, m, fp, b, xi, n0, init.model=init.model, method="update.Q")

    fp <- BHLK.filter(SS$e, Xik, Q)

    fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize

#    plot(ts(fp), plot.type="single")

    # Print interim outputs

    cat("Iteration : ", iter, "\n")
    cat("b  : ", b, "\n")
#    cat("xi : ", xi, "\n")
    cat("Q  : ", Q, "\n")


    if(iter>1)
    {
        cat("Relative tol. b  : ", sum(b-b.old)/sum(b.old^2), "\n")
#        cat("Relative tol. xi : ", sum(xi-xi.old)/sum(xi.old^2), "\n")
        cat("Relative tol. Q  : ", sum(Q-Q.old)/sum(Q.old^2), "\n")
        cat("\n")

        cat("Relative fn. b   : ", (max.b$value-max.b.fn)/(max.b.fn), "\n")
#        cat("Relative fn. xi  : ", (max.xi$value-max.xi.fn)/(max.xi.fn), "\n")
        cat("Relative fn. Q   : ", (max.Q$value-max.Q.fn)/(max.Q.fn), "\n")


        cat("Function value b  :", max.b$value, "\n")
#        cat("Function value xi :", max.xi$value, "\n")
        cat("Function value Q  :", max.Q$value, "\n")

        # Reset function values
        max.b.fn <- max.b$value
#        max.xi.fn <- max.xi$value
        max.Q.fn <- max.Q$value

    }
}

    # Define the output object
    output <- list(b=b,
                   xi=xi,
                   Q=Q,
                   fp=fp,
                   m=m, p=p, h=h,
                   init.model=init.model,
                   n0=n0,
                   n0cum=n0cum,
                   abar=abar,
                   bbar=bbar,
                   alpha.prior=alpha.prior)

    class(output) <- c("MSBSVARsetup")

    return(output)
}

