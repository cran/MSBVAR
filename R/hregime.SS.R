# H-regime SS regressions by equation -- can use this to get
# regression normal equations.
hregime.SS <- function(h, m, fp, b, xi, n0, init.model, method=NULL)
{

    # Subset init.model variables
    Y <- init.model$Y[(m+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1):nrow(init.model$X),]
    p <- dim(init.model$ar.coefs)[3]

    n0cum <- c(0, cumsum(n0))
    # split xi and b by regime
    xi <- matrix(xi, m, h)
    b <- matrix(b, ncol=h)

    # normalize xi
    xi <- norm.xi(xi)

    # Find obs. for each state
    fp <- round(fp)

#    fp <- fp[,c(rev(rank(colSums(fp), ties.method="first")))]   #normalize

    # Storage
    tmp <- vector(mode="list", length=h)
    bjk.cov <- vector(mode="list", length=h)

    for(i in 1:h)
    { bjk.cov[[i]] <- vector(mode="list", length=m)
    }

    A0k <- array(0, c(m,m,h))
    F <- array(0, c(m*p+1, m, h))
    e <- array(0, c(nrow(fp), m, h))

    Sbar <- rbind(diag(m), matrix(0, (m*(p-1)+ 1), m))

    # get dfs for each regime
#    df <- colSums(fp) - m*p + 1

    # Loops to compute
    # 1) Sums of squares for X and Y
    # 2) Number of observations in each state.
    # 3) A0k matrices
    # 4) F(k) matrices
    # 5) Residuals e(k)

    for(i in 1:h)
    {
        A0k[,,i] <- b2a(b[,i], init.model$Ui)

        for (j in 1:m)
        {
            ytmp <- matrix(Y[which(fp[,i]==1),], ncol=m)
            xtmp <- matrix(X[which(fp[,i]==1),], ncol=ncol(X))

            Syy <- crossprod(ytmp)*xi[j,i]
            Sxy <- crossprod(xtmp, ytmp)*xi[j,i]
            Sxx <- crossprod(xtmp)*xi[j,i]

            tmp[[i]] <- list(nrow(ytmp), Syy, Sxy, Sxx, ytmp, xtmp)

            fjk <- Sxx + init.model$Hpinv.tilde[,,j]
            XY.post <- Sxy%*%init.model$Ui[[j]] + init.model$Hpinv.tilde[,,j]%*%init.model$Pi.tilde[[j]]
            P.post <- solve(fjk, XY.post)
            btmp <- b[(n0cum[j]+1):(n0cum[(j+1)]),i]
            F[,j,i] <- P.post%*%btmp


            if (method=="update.b") {
            # posterior covariance of bjk -- structural coef. for
            # equation j in state k
##             bjk.cov[[i]][[j]] <- solve(t(init.model$Ui[[j]])%*%Syy%*%init.model$Ui[[j]] +
##                                        init.model$H0inv.tilde[[j]] +
##                                        t(init.model$Pi.tilde[[j]])%*%init.model$Hpinv.tilde[,,j]%*%init.model$Pi.tilde[[j]] -
##                                        t(XY.post)%*%P.post)

                bjk.cov[[i]][[j]] <- solve(init.model$H0inv.tilde[[j]] +
                                           crossprod(init.model$Ui[[j]], (Syy + crossprod(Sbar,Sxy) + crossprod(Sxy, Sbar) +
                                                                          crossprod(Sbar, Sxx)%*%Sbar))%*%init.model$Ui[[j]])


            }

        }

        if (method=="update.Q") { e[,,i] <- Y%*%A0k[,,i] - X%*%F[,,i]
                              }

    }
    return(list(moments=tmp, A0k=A0k, F=F, e=e, bjk.cov=bjk.cov))
}

