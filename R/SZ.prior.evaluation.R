"SZ.prior.evaluation" <- function(Y, p, lambda0, lambda1, lambda3,
  lambda4, lambda5, mu5, mu6, z=NULL, nu=ncol(Y)+1, qm, prior=0,
  nstep, y.future)
  { combos <- length(lambda0)*length(lambda1)*length(lambda3)*length(lambda4)
    results <- matrix(0, combos ,11)
    h <- 0

    for(i in 1:length(lambda0))
      { for(j in 1:length(lambda1))
        { for(k in 1:length(lambda3))
          { for(l in 1:length(lambda4))
            { fit <- szbvar(Y, p, z,
                            lambda0[i],
                            lambda1[j],
                            lambda3[k],
                            lambda4[l],
                            lambda5=lambda5,
                            mu5=mu5,
                            mu6=mu6,
                            nu=nu,
                            qm=qm,
                            prior=prior,
                            posterior.fit=T)

              forecast <- forecast.VAR(fit, nstep)
              eval.forecasts <-
                cf.forecasts(forecast[(nrow(forecast)-nstep+1):nrow(forecast),], y.future)

              # Now gather up all the parameters and results
              tmp <- c(lambda0[i],lambda1[j],lambda3[k],lambda4[l],
                       lambda5, mu5, mu6, eval.forecasts[1],
                       eval.forecasts[2], fit$marg.llf[1],
                       fit$marg.post[1])

              h <- h+1
              results[h,] <- tmp

            }}}}

    colnames(results) <- c("lambda0", "lambda1", "lambda3", "lambda4",
                           "lambda5", "mu5", "mu6", "RMSE", "MAE", "MargLLF",
                           "MargPosterior")

    return(results)
  }
