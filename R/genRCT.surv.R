#' Generalizing the average treatment effect (ATE) from the trial using observational studies (OS) for survival outcomes
#'
#' Provides estimation and inference for the average treatment effect (ATE), generalizing from a randomized controlled trial (RCT)
#' to a target population leveraging observational studies. The method of sieves can be used for the estimation of sampling score and outcome models.
#'
#' @param Y.trial Observed outcome from a trial; vector of size \code{n} (trial sample size).
#' @param d.trial The event indicator, normally 1 = event, 0 = censored; vector of size \code{n}.
#' @param A.trial Treatment received from a trial; vector of size \code{n}.
#' @param X.trial Matrix of \code{p} baseline covariates from a trial; dimension \code{n} by \code{p}.
#' @param X.rwe Matrix of \code{p} baseline covariates from OS; dimension \code{m} by \code{p}.
#' @param tau | A vector of truncation time for defining restricted mean survival time; e.g., seq(10, 50, by = 10)
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{100}.
#' @param conf.level The level of bootstrap confidence interval; Default is \code{0.05}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose A logical value indicating whether intermediate progress messages should be printed.
#' Default is \code{TRUE}.
#' @return Return a list containing:\tabular{ll}{
#'    \code{rmst} \tab\tab  A list of estimated RMSTs with bootstrap SE and confidence interval. \cr
#'    \tab \cr
#'    \code{surv} \tab\tab  A list of estimated treatment-specific survival function. \cr
#' }
#'
#' @examples
#' fit <- genRCT.surv(Y.trial = Y.trial, d.trial = d.trial, A.trial = A.trial, X.trial = X.trial, X.rwe = X.rwe,
#' tau = tau. n.boot = 100, conf.level = 0.05, seed = 123, verbose = TRUE)
#' @export

'genRCT.surv' <- function(Y.trial, d.trial, A.trial, X.trial, X.rwe, tau, eta.vec = c(0,0.001,0.005),
                          n.boot = 100, conf.level = 0.05, seed = NULL, verbose = TRUE) {

  t.start <- Sys.time()

  # Change data.frame to matrix
  Y.trial <- c(as.matrix(Y.trial))
  A.trial <- c(as.matrix(A.trial))
  X.trial <- as.matrix(X.trial)
  X.rwe <- as.matrix(X.rwe)

  # Fit estimators
  cat(" Fitting estimators.. \n")
  fit <- genRCT.surv.estimators(Y.trial = Y.trial, d.trial = d.trial, A.trial = A.trial, X.trial = X.trial, X.rwe = X.rwe, O.sel1 = NULL, O.sel0 = NULL,
                                 C.sel1 = NULL, C.sel0 = NULL, A.sel = NULL, eta = NULL, eta.vec = eta.vec, seed = seed)

  time1 <- sort(Y.trial[A.trial == 1 & d.trial == 1])
  time0 <- sort(Y.trial[A.trial == 0 & d.trial == 1])
  time1.idx <- which(fit$surv1$time %in% time1)
  time0.idx <- which(fit$surv0$time %in% time0)
  t0 <- c(0, rep(1, 14))
  surv1 <- rbind(t0, fit$surv1[time1.idx, ])
  surv0 <- rbind(t0, fit$surv0[time0.idx, ])

  fit.rmst <- rmst_to_save(fit.obj = fit, tau = tau)


  n1 <- nrow(X.trial)
  n2 <- nrow(X.rwe)
  survB.1 <- array(NA, dim = c(nrow(surv1), ncol(surv1), n.boot))
  survB.0 <- array(NA, dim = c(nrow(surv0), ncol(surv0), n.boot))
  boot.rmst <- array(NA, dim = c(nrow(fit.rmst), ncol(fit.rmst), n.boot))

  if(is.null(seed)) seed <- 0
  cat(" Bootstrapping.. \n")
  if(verbose == TRUE) {
    pb <- txtProgressBar(min = 0, max = n.boot, style = 3)
  }
  for (i in 1:n.boot){
    #samp1b <- sample(n1, replace = TRUE)
     #sampA1 <- sample(which(A.trial == 1), replace = TRUE)
     #sampA0 <- sample(which(A.trial == 0), replace = TRUE)
     #samp1b <- c(sampA1, sampA0)
     samp1b <- sample(n1, replace = TRUE)
     samp2b <- sample(n2, replace = TRUE)

    #sampA1 <- balance.boots(x = X.trial, A.trial, trt = 1, rwe = FALSE)
    #sampA0 <- balance.boots(x = X.trial, A.trial, trt = 0, rwe = FALSE)
    #samp1b <- c(sampA1, sampA0)
    #samp2b <- balance.boots(x = X.rwe, A.trial, trt = 0, rwe = TRUE)

    X1b <- X.trial[samp1b, ]
    d1b <- d.trial[samp1b]
    A1b <- A.trial[samp1b]
    Y1b <- Y.trial[samp1b]
    X2b <- X.rwe[samp2b, ]
    fit.boot <- try(genRCT.surv.estimators(Y.trial = Y1b, d.trial = d1b, A.trial = A1b, X.trial = X1b, X.rwe = X2b,
                                            O.sel1 = fit$hyper$O.sel1, O.sel0 = fit$hyper$O.sel0, C.sel1 = fit$hyper$C.sel1, C.sel0 = fit$hyper$C.sel0,
                                            A.sel = fit$hyper$A.sel, eta = fit$hyper$eta, eta.vec = eta.vec, seed = seed + i), silent = TRUE)

    if(inherits(fit.boot, "try-error")) next;

    time1.all <- union(time1, fit.boot$surv1$time)
    time0.all <- union(time0, fit.boot$surv0$time)

    s1.B <- rbind(t0, data.frame(time = time1.all) %>% left_join(fit.boot$surv1, by = c("time"))) %>% fill(surv.nv:denom.dacw)
    s0.B <- rbind(t0, data.frame(time = time0.all) %>% left_join(fit.boot$surv0, by = c("time"))) %>% fill(surv.nv:denom.dacw)

    survB.1[, , i] <- as.matrix(s1.B %>% filter(time %in% c(0, time1)))
    survB.0[, , i] <- as.matrix(s0.B %>% filter(time %in% c(0, time0)))
    boot.rmst[, , i] <- rmst_to_save(fit.obj = fit.boot, tau = tau)

    if (verbose == TRUE) {
      Sys.sleep(0.025)
      setTxtProgressBar(pb, i)
    }
  } # end of for(i in 1:n.boot)
  if (verbose == TRUE) close(pb)

  se1 <- apply(survB.1, c(1, 2), sd, na.rm = TRUE)
  se0 <- apply(survB.0, c(1, 2), sd, na.rm = TRUE)
  lower1 <- apply(survB.1, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  lower0 <- apply(survB.0, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  upper1 <- apply(survB.1, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  upper0 <- apply(survB.0, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)

  colnames(se1) <- colnames(se0) <- colnames(surv1)
  colnames(lower1) <- colnames(lower0) <- colnames(upper1) <- colnames(upper0) <- colnames(surv1)

  se.rmst <- apply(boot.rmst, c(1, 2), sd, na.rm = TRUE)
  lower.rmst <- apply(boot.rmst, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  upper.rmst <- apply(boot.rmst, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)

  colnames(se.rmst) <- colnames(lower.rmst) <- colnames(lower.rmst) <- colnames(fit.rmst)
  rownames(se.rmst) <- rownames(lower.rmst) <- rownames(upper.rmst) <- rownames(fit.rmst)

  surv <- list(time1 = surv1$time, surv1 = surv1[, -1], se1 = se1[,-1], lower1 = lower1[,-1], upper1 = upper1[,-1],
               time0 = surv0$time, surv0 = surv0[, -1], se0 = se1[,-1], lower0 = lower0[,-1], upper0 = upper0[,-1])
  rmst <- list(est = fit.rmst, se = se.rmst, lower = lower.rmst, upper = upper.rmst)

  out <- c()
  out$surv <- surv
  out$rmst <- rmst
  out$hyper <- fit$hyper

  t.end <- Sys.time()
  t.diff <- difftime(t.end, t.start)
  cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')

  return(out)

}



