# plotregimeid() : clustering and plotting function for msbvar permuted
#                   sample output
#
# Iliyan Iliev and Patrick T. Brandt
#
# 20101115 : Initial version
# 20101124 : Added more plots and first logical switches for types of
#            plots for posterior identification
# 20101130 : Added labelings for plots; begin adding logical
#            conditions for the plot types.
# 20110118 : Updated logical conditions for plots and added section
#            for plotting versus values of Q.

plotregimeid <- function(m1, m1p,
                         type=c("all", "intercepts", "Sigma", "Q"),
                         ask = TRUE, ...)
{

    devAskNewPage(ask)
    # Get the constants we need for setting up the matrices for the
    # clustering steps.

    m <- m1$m                     # Number of equations
    p <- m1$p                     # Lag length
    h <- m1$h                     # Number of regimes
    N2 <- length(m1p$ss.sample)   # Number of Gibbs draws

    mpplus1 <- m*p + 1  # Number of coefficients in one equation of
                        # VAR(p)

    nc <- m*(mpplus1)   # Number of total coefs in one regime

    ii <- seq(mpplus1, by=mpplus1, length=m)  # Intercept indices
    eqnnames <- colnames(m1$init.model$y)

    if(type=="all" || type=="intercepts")
    {
    # Now stack the regression coefficients for clustering
        Beta <- matrix(t(m1p$Beta.sample), byrow=TRUE, ncol=nc)

    # Now do the clustering for the regimes.
    # Probably should let the users choose the number of starting
    # points for the clustering centers = nstart, later.

        intercepts <- as.data.frame(Beta[,ii])
    # Label things so the plots make sense
        colnames(intercepts) <- eqnnames

        cl.int <- kmeans(intercepts, center=h, nstart=10)


    # Plot the intercepts based on the clustering
        pairs(intercepts, pch=".", col=cl.int$cluster)
        title("Intercepts pairs by regime", line=3)

    # Note that the lattice / coda plotting function (the first two),
    # need to wrapped in a print to be displayed.

	form <- as.formula(paste("~",
                                 paste(lapply(names(intercepts), as.name),
                                       collapse = "+")))

        print(densityplot(form, data=intercepts, outer=TRUE,
                    groups=cl.int$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Intercept densities by regime"))

        form <- eval(parse(text = paste(paste(lapply(names(intercepts), as.name),
                           collapse = "+"), "~ idx")))
        idx <- 1:nrow(intercepts)

        print(xyplot(form, data=intercepts,
                     groups=cl.int$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Intercept traceplots by regime"))
    }

    if(type=="all" || type=="Sigma")
    {

    # Now replicate this for the variance elements in
    # m1p$Sigma.sample.  Note that the unique elements of the Sigma(h)
    # matrices are in this and make up a different dimension than the
    # regression effects.
        nvar <- (m*(m+1)*0.5) # Number of variance terms
        Sigmaout <- matrix(t(m1p$Sigma), byrow=TRUE, ncol=nvar)

    # Get variance indices
        tmp <- m:2
        tmp <- rep(c(1,tmp), h)
        ssidx <- tmp
        for(i in 2:(m*h)) ssidx[i] <- ssidx[i-1] + tmp[i]

        Sigmaout <- as.data.frame(Sigmaout[,ssidx[1:m]])
        colnames(Sigmaout) <- eqnnames

        cl.sigma <- kmeans(Sigmaout, center=h, nstart=10)

    # Plots
        pairs(Sigmaout, pch=".", col=cl.sigma$cluster, ...)
        title("Variances pairs plot by regime", line=3)

        form <- as.formula(paste("~",
                          paste(lapply(names(Sigmaout), as.name),
                                collapse = "+")))

        print(densityplot(form, data=Sigmaout, outer=TRUE,
                    groups=cl.sigma$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Variance densities by regime"), ...)

        form <- eval(parse(text = paste(paste(lapply(names(Sigmaout), as.name),
                           collapse = "+"), "~ idx")))
        idx <- 1:nrow(Sigmaout)

        print(xyplot(form, data=Sigmaout,
                     groups=cl.sigma$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Variance traceplots by regime"), ...)
    }

    if(type=="all" || type=="Q")
    {
    # Do the same as the above for the elements of m1p$Q.
    # These do not need to be stacked like the other elements (why?)

        Q <- as.data.frame(m1p$Q)

    # Cluster

        cl.Q <- kmeans(Q, center=h, nstart=10)

        idx <- cbind(rep(1:h, each=h), rep(1:h, times=h))
        Qnames <- paste("Q_", idx[,1], idx[,2], sep="")
        colnames(Q) <- Qnames

    # Plots

        pairs(Q, pch=".", col=cl.Q$cluster, ...)
        title("Transitions pairs plot by regime", line=3)

	form <- as.formula(paste("~",
                          paste(lapply(names(Q), as.name),
                                collapse = "+")))

        print(densityplot(form, data=Q, outer=TRUE,
                    groups=cl.Q$cluster, xlab=NULL,
                    default.scales=list(relation="free"),
                    plot.points="rug",
                    main="Transition densities by regime"), ...)

        form <- eval(parse(text = paste(paste(lapply(names(Q), as.name),
                           collapse = "+"), "~ idx")))
        idx <- 1:nrow(Q)

        print(xyplot(form, data=Q,
                     groups=cl.Q$cluster, type="l", ylab=NULL,
                     default.scales=list(relation="free"),
                     main="Transition traceplots by regime"), ...)


    # Combine all of the above into one set of plots / analyses across
    # intercepts -- like SFS 2001, Figure 8

    # Plots based on the posterior clustering of Q for the intercepts
       # Now stack the regression coefficients for clustering
        Beta <- matrix(t(m1p$Beta.sample), byrow=TRUE, ncol=nc)

    # Now do the clustering for the regimes.
    # Probably should let the users choose the number of starting
    # points for the clustering centers = nstart, later.

        intercepts <- as.data.frame(Beta[,ii])
    # Label things so the plots make sense
        colnames(intercepts) <- eqnnames

    # Get diagonal of Q
        qdiag <- diag(matrix(1:h^2, h, h))

        par(mfrow=c(2, round(m/2)), omi=c(0.5, 0.75, 0.75, 0.25))
        for(i in 1:m)
        {
            plot(intercepts[,i], matrix(unlist(Q[, qdiag]), ncol=1), pch=".",
                 col=cl.Q$cluster, xlab=names(intercepts)[i],
                 ylab="Transition Probability Regimes")
        }
        title("Intercepts by transition probability regimes",
              outer=TRUE, line=1)

    }
    devAskNewPage(FALSE)
}
