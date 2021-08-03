#' Estimate the ATE using specified estimators
#'
#' ATE estimated with Naive, IPSW, AIPSW, CW, ACW-t, ACW-b.
#'
#' @param Y.trial Observed outcome from a trial; vector of size \code{n} (the trial sample size).
#' @param A.trial Treatment received from a trial; vector of size \code{n}.
#' @param X.trial Matrix of \code{p} baseline covariates from a trial; dimension \code{n} by \code{p}.
#' @param Y.rwe Observed outcome from OS; if obtained, vector of size \code{m} (OS sample size);
#' otherwise, set \code{Y.rwe = NULL}.
#' @param A.rwe Treatment received from OS; if obtained, vector of size \code{m}; otherwise, set \code{Y.rwe = NULL}.
#' @param X.rwe Matrix of \code{p} baseline covariates from OSl; dimension \code{m} by \code{p}.
#' @param family The type of outcome; \code{"gaussian"} for continuous outcome or \code{"binomial"} for binary outcome.
#' Default is \code{"gaussian"}.
#' @param estimators A vector of one or multiple methods to estimate the ATE. Allowed values are
#' \code{'Naive'}, \code{'IPSW'}, \code{'AIPSW'}, \code{'CW'}, \code{'ACW-t'}, \code{'ACW-b'}.
#' The \code{'ACW-b'} is allowed only when both \code{Y.rwe} and \code{A.rwe} are obtained.
#' Default specifies all 6 methods.
#' @param sieve A logical value indicating whether the method of sieves are used for estimating sampling score and outcome models.
#' Used only if \code{estimators = 'AIPSW'} or \code{'ACW-t'} or \code{'ACW-b'}. Default is \code{TRUE}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param osel1.t A vector indicating selected outcome model covariates for the \code{'ACW-t'} estimator in the group trt = 1.
#' This is only relevant when \code{sieve = TRUE} and \code{inference = TRUE}. Otherwise, set \code{osel1.t = NULL}.
#' @param osel0.t A vector indicating selected outcome model covariates for the \code{'ACW-t'} estimator in the group trt = 0.
#' This is only relevant when \code{sieve = TRUE} and \code{inference = TRUE}. Otherwise, set \code{osel0.t = NULL}.
#' @param osel1.b A vector indicating selected outcome model covariates for the \code{'ACW-b'} estimator in the group trt = 1.
#' This is only relevant when \code{sieve = TRUE} and \code{inference = TRUE}. Otherwise, set \code{osel1.b = NULL}.
#' @param osel0.b A vector indicating selected outcome model covariates for the \code{'ACW-b'} estimator in the group trt = 0.
#' This is only relevant when \code{sieve = TRUE} and \code{inference = TRUE}. Otherwise, set \code{osel0.b = NULL}.
#' @param osel.ipsw A vector indicating selected sampling scole model covariates for the \code{'AIPSW'} estimator.
#' This is only relevant when \code{sieve = TRUE} and \code{inference = TRUE}. Otherwise, set \code{osel.ipsw = NULL}.
#' @return Return a list containing:\tabular{ll}{
#'    \code{ate} \tab  A vector of estimated ATEs using the selected methods. \cr
#'    \tab \cr
#'    \code{hypers} \tab A list of vectors indicating the selected covariates. \cr
#' }
#'

genRCT.estimators <- function(Y.trial, A.trial, X.trial, Y.rwe, A.rwe , X.rwe, family = "gaussian", estimators, sieve = TRUE,  seed = NULL,
                              osel1.t = NULL, osel0.t = NULL, osel1.b = NULL, osel0.b = NULL, osel.ipsw = NULL) {
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Names
  ate <- matrix(NA, nrow = 1, ncol = length(estimators))
  colnames(ate) <- c(paste0(estimators))
  ate <- data.frame(ate)

  dat.trial <- data.table(A = A.trial, X = X.trial)
  dat.trial$Y <- Y.trial

  if (!is.null(Y.rwe) & !is.null(A.rwe)) {
    dat.rwe <- data.table(A = A.rwe, X = X.rwe)
    dat.rwe$Y <- Y.rwe
  } else {
    dat.rwe <- data.table(X = X.rwe)
  }

  n.trial <- nrow(X.trial)
  n.rwe <- nrow(X.rwe)
  p.trial <- ncol(X.trial)
  p.rwe <- ncol(X.rwe)
  X.all <- rbind(X.trial, X.rwe)

  p1 <- sum(dat.trial$A == 1) / n.trial
  p0 <- 1 - p1
  delta <- c(rep(1, n.trial), rep(0, n.rwe))

  # Fit outcome regression

  H.trial <- cbind(X.trial, X.trial * A.trial)
  h1.trial <- cbind(1, X.trial, X.trial * 1)
  h0.trial <- cbind(1, X.trial, X.trial * 0)
  h1.rwe <- cbind(1, X.rwe, X.rwe * 1)
  h0.rwe <- cbind(1, X.rwe, X.rwe * 0)

  if (sieve == TRUE) {
    gvec.trial <- model.matrix(~ poly(X.trial, degree = 2, raw = TRUE))
    gvec.trial <- gvec.trial[, which(!duplicated(t(gvec.trial))), drop = FALSE]
    gvec.rwe <- model.matrix(~ poly(X.rwe, degree = 2, raw = TRUE))
    gvec.rwe <- gvec.rwe[, which(!duplicated(t(gvec.rwe))), drop = FALSE]
  }

  if (any(c("AIPSW", "ACW-t") %in% estimators)) {

    if (sieve == FALSE) {

      if(family == 'gaussian') {
        fit.t <- lm(dat.trial$Y ~ H.trial)
        beta.t <- as.numeric(coef(fit.t))
        dat.trial$Y1.t <- h1.trial %*% beta.t
        dat.trial$Y0.t <- h0.trial %*% beta.t
        dat.rwe$Y1.t <- h1.rwe %*% beta.t
        dat.rwe$Y0.t <- h0.rwe %*% beta.t
      } else if (family == 'binomial') {
        fit.t <- glm(dat.trial$Y ~ H.trial, family = binomial("logit"))
        beta.t <- as.numeric(coef(fit.t))
        dat.trial$Y1.t <- expit(h1.trial %*% beta.t)
        dat.trial$Y0.t <-  expit(h0.trial %*% beta.t)
        dat.rwe$Y1.t <-  expit(h1.rwe %*% beta.t)
        dat.rwe$Y0.t <-  expit(h0.rwe %*% beta.t)
      }
    } else {

      designX <- gvec.trial[, 2:ncol(gvec.trial)]
      indA1.trial <- which(dat.trial$A == 1)
      indA0.trial <- which(dat.trial$A == 0)

      if (is.null(osel1.t) | is.null(osel0.t)) {
        for(ns in 1:5) {
          fit1.t.sv <- cv.ncvreg(designX[indA1.trial, ], dat.trial$Y[indA1.trial], penalty = 'SCAD', family = family)#, warn=FALSE)
          beta1.t.sv <- as.numeric(coef(fit1.t.sv, s = 'lambda.min'))
          osel1.t <- union(osel1.t, which(beta1.t.sv != 0))
          fit0.t.sv <- cv.ncvreg(designX[indA0.trial, ], dat.trial$Y[indA0.trial], penalty = 'SCAD', family = family)#, warn=FALSE)
          beta0.t.sv <- as.numeric(coef(fit0.t.sv, s = 'lambda.min'))
          osel0.t <- union(osel0.t, which(beta0.t.sv != 0))
        }
      }

      #Refit the selected basis
      if(family == 'gaussian'){
        fit1.t.sv <- lm(dat.trial$Y[indA1.trial] ~ cbind(1, designX)[indA1.trial, osel1.t] - 1)
        beta1.t.sv <- as.numeric(coef(fit1.t.sv))
        dat.trial$Y1.t.sv <- gvec.trial[, osel1.t, drop = FALSE] %*% beta1.t.sv
        dat.rwe$Y1.t.sv <- gvec.rwe[, osel1.t, drop = FALSE] %*% beta1.t.sv
        fit0.t.sv <- lm(dat.trial$Y[indA0.trial] ~ cbind(1, designX)[indA0.trial, osel0.t] - 1)
        beta0.t.sv <- as.numeric(coef(fit0.t.sv))
        dat.trial$Y0.t.sv <- gvec.trial[, osel0.t, drop = FALSE] %*% beta0.t.sv
        dat.rwe$Y0.t.sv <- gvec.rwe[, osel0.t, drop = FALSE] %*% beta0.t.sv
      } else if (family == 'binomial') {
        fit1.t.sv <- glm(dat.trial$Y[indA1.trial] ~ cbind(1, designX)[indA1.trial, osel1.t] - 1, family = binomial("logit"))
        beta1.t.sv <- as.numeric(coef(fit1.t.sv))
        dat.trial$Y1.t.sv <- expit(gvec.trial[, osel1.t, drop = FALSE] %*% beta1.t.sv)
        dat.rwe$Y1.t.sv <- expit(gvec.rwe[, osel1.t, drop = FALSE] %*% beta1.t.sv)
        fit0.t.sv <- glm(dat.trial$Y[indA0.trial] ~ cbind(1, designX)[indA0.trial, osel0.t] - 1, family = binomial("logit"))
        beta0.t.sv <- as.numeric(coef(fit0.t.sv))
        dat.trial$Y0.t.sv <- expit(gvec.trial[, osel0.t, drop = FALSE] %*% beta0.t.sv)
        dat.rwe$Y0.t.sv <- expit(gvec.rwe[, osel0.t, drop = FALSE] %*% beta0.t.sv)
      }
    } # end of if (sieve == FALSE) {} else {}
  } # end of if (any("AIPSW", "ACW-t") %in% estimators)

  if ("ACW-b" %in% estimators) {

    H.rwe <- cbind(X.rwe, X.rwe * A.rwe)
    if (sieve == FALSE) {
      if(family == 'gaussian') {
        fit.b <- lm(dat.rwe$Y ~ H.rwe)
        beta.b <- as.numeric(coef(fit.b))
        dat.trial$Y1.b <- h1.trial %*% beta.b
        dat.trial$Y0.b <- h0.trial %*% beta.b
        dat.rwe$Y1.b <- h1.rwe %*% beta.b
        dat.rwe$Y0.b <- h0.rwe %*% beta.b
      } else if (family == 'binomial') {
        fit.b <- glm(dat.rwe$Y ~ H.rwe, family = binomial("logit"))
        beta.b <- as.numeric(coef(fit.b))
        dat.trial$Y1.b <- expit(h1.trial %*% beta.b)
        dat.trial$Y0.b <- expit(h0.trial %*% beta.b)
        dat.rwe$Y1.b <- expit(h1.rwe %*% beta.b)
        dat.rwe$Y0.b <-  expit(h0.rwe %*% beta.b)
      }

    } else {

      designX <- gvec.rwe[, 2:ncol(gvec.rwe)]
      indA1.rwe <- which(dat.rwe$A == 1)
      indA0.rwe <- which(dat.rwe$A == 0)

      if(is.null(osel1.b) | is.null(osel0.b)){
        for(ns in 1:5){
          fit1.b.sv <- cv.ncvreg(designX[indA1.rwe, ], dat.rwe$Y[indA1.rwe], penalty = "SCAD", family = family)#, warn=FALSE)
          beta1.b.sv <- as.numeric(coef(fit1.b.sv, s = "lambda.min"))
          osel1.b <- union(osel1.b, which(beta1.b.sv != 0))
          fit0.b.sv <- cv.ncvreg(designX[indA0.rwe, ], dat.rwe$Y[indA0.rwe], penalty = "SCAD", family = family)#, warn=FALSE)
          beta0.b.sv <- as.numeric(coef(fit0.b.sv, s = "lambda.min"))
          osel0.b <- union(osel0.b, which(beta0.b.sv != 0))
        }
      }

      #Refit the selected basis
      if (family == 'gaussian'){
        fit1.b.sv <- lm(dat.rwe$Y[indA1.rwe] ~ cbind(1, designX)[indA1.rwe, osel1.b] - 1)
        beta1.b.sv <- as.numeric(coef(fit1.b.sv))
        dat.trial$Y1.b.sv <- gvec.trial[, osel1.b, drop = FALSE] %*% beta1.b.sv
        dat.rwe$Y1.b.sv <- gvec.rwe[, osel1.b, drop = FALSE] %*% beta1.b.sv
        fit.0.sv <- lm(dat.rwe$Y[indA0.rwe] ~ cbind(1, designX)[indA0.rwe, osel0.b] - 1)
        beta0.b.sv <- as.numeric(coef(fit.0.sv))
        dat.trial$Y0.b.sv <- gvec.trial[, osel0.b, drop = FALSE] %*% beta0.b.sv
        dat.rwe$Y0.b.sv <- gvec.rwe[, osel0.b, drop = FALSE] %*% beta0.b.sv
      } else if (family == 'binomial') {
        fit1.b.sv <- glm(dat.rwe$Y[indA1.rwe] ~ cbind(1, designX)[indA1.rwe, osel1.b] - 1, family = binomial("logit"))
        beta1.b.sv <- as.numeric(coef(fit1.b.sv))
        dat.trial$Y1.b.sv <- expit(gvec.trial[, osel1.b, drop = FALSE] %*% beta1.b.sv)
        dat.rwe$Y1.b.sv <- expit(gvec.rwe[, osel1.b, drop = FALSE] %*% beta1.b.sv)
        fit.0.sv <- glm(dat.rwe$Y[indA0.rwe] ~ cbind(1, designX)[indA0.rwe, osel0.b] - 1, family = binomial("logit"))
        beta0.b.sv <- as.numeric(coef(fit.0.sv))
        dat.trial$Y0.b.sv <- expit(gvec.trial[, osel0.b, drop = FALSE] %*% beta0.b.sv)
        dat.rwe$Y0.b.sv <- expit(gvec.rwe[, osel0.b, drop = FALSE] %*% beta0.b.sv)
     }
   } # end of if (sieve == FALSE) {} else {}
 } # end of if ("ACW-b" %in% estimators)


 # Estimate calibration weights
  if ("CW" %in% estimators | (any(c("ACW-t", "ACW-b") %in% estimators & sieve == FALSE))) {

    moment.bar <- colMeans(X.rwe)
    moms.t <- cbind(X.trial)

    cnt <- 0
    while(cnt <= 50){
      cnt <- cnt + 1
      lam.hat <- searchZeros(matrix(rnorm(length(moment.bar) * 20, 0, 0.25), nrow = 20), lamFun, moments = X.t, moments.bar = moment.bar)$x[1,]
      if(!is.null(lam.hat)) break
    }

    if(is.null(lam.hat)){
      warning('No lam.hat solutions')
      lam.hat <- rep(NA, p.trial)
    }
    q.score <- exp(X.trial %*% lam.hat) / sum(exp(X.trial %*% lam.hat))
    dat.trial$q <- q.score
  } # end of if ( "CW" %in% estimators | (any(c("ACW-t", "ACW-b") %in% estimators & sieve == FALSE)))

  if (c("ACW-t") %in% estimators & sieve == TRUE) {
    osel.t <- union(osel1.t, osel0.t)
    if(length(osel.t) > 1){
      q.score <- backward.sv(osel1 = osel1.t, osel0 = osel0.t, beta1.sv = beta1.t.sv, beta0.sv = beta0.t.sv, gvec.trial = gvec.trial, gvec.rwe = gvec.rwe)
    } else {
      q.score <- 1 / n.trial
    }
    dat.trial$q.t.sv <- q.score
  }

  if (c("ACW-b") %in% estimators & sieve == TRUE) {
    osel.b <- union(osel1.b, osel0.b)
    if(length(osel.b) > 1){
      q.score <- backward.sv(osel1 = osel1.b, osel0 = osel0.b, beta1.sv = beta1.b.sv, beta0.sv = beta0.b.sv, gvec.trial = gvec.trial, gvec.rwe = gvec.rwe)
    } else {
      q.score <- 1 / n.trial
    }
    dat.trial$q.b.sv <- q.score
  }


  # Estimate (A)IPSW weight
  if (any(c("IPSW", "AIPSW") %in% estimators)) {
    fit.ipsw <- glm(delta ~ X.all, family = binomial("logit"))
    phi.hat <- fit.ipsw$coefficients
    s.score <- as.numeric(expit(cbind(1, X.trial) %*% phi.hat))
    s.score <- pmin(s.score, 0.99)
    s.score <- pmax(s.score, 0.01)
    dat.trial$ipsw.w <- 1 / s.score
  }

  ##### Fit estimators #####

  if ("Naive" %in% estimators) {
    tau.naive <- mean(dat.trial[dat.trial$A == 1, ]$Y) - mean(dat.trial[dat.trial$A == 0, ]$Y)
    ate$Naive <- tau.naive
  }

  if ("IPSW" %in% estimators) {
    tau.ipsw <- with(dat.trial, sum((A * Y * ipsw.w) / sum(A * ipsw.w) - (((1 - A) * Y * ipsw.w)) / sum((1 - A) * ipsw.w)))
    ate$IPSW <- tau.ipsw
  }

  if ("AIPSW" %in% estimators) {
    if (sieve == FALSE) {
      tau.aipsw <- with(dat.trial, sum(ipsw.w * ((A * (Y - Y1.t) / sum(A * ipsw.w) - ((1 - A) * (Y - Y0.t)) / sum((1 - A) * ipsw.w))))) + with(dat.rwe, mean(Y1.t - Y0.t))
    } else {
      gvec.all <- model.matrix( ~ poly(X.all,degree = 2,raw = TRUE))
      gvec.all <- gvec.all[, which(!duplicated(t(gvec.all))), drop = FALSE]

      designX <- gvec.all[, 2:ncol(gvec.all)]
      if (is.null(osel.ipsw)) {
        for(ns in 1:5){
          fit.ipsw.sv <- cv.ncvreg(designX, delta, penalty = "SCAD", family = "binomial", warn = TRUE)
          phi.hat.sv <- as.numeric(coef(fit.ipsw.sv, s = "lambda.min"))
          osel.ipsw <- union(osel.ipsw, which(phi.hat.sv != 0))
        }
      }
      # Refit the selected basis
      fit.ipsw.sv <- glm(delta ~ cbind(1, designX)[, osel.ipsw] - 1, family = binomial("logit"))
      phi.hat.sv <- fit.ipsw.sv$coefficients
      s.score.sv <- as.numeric(expit(gvec.trial[, osel.ipsw, drop = FALSE] %*% phi.hat.sv))
      s.score.sv <- pmin(s.score.sv, 0.99)
      s.score.sv <- pmax(s.score.sv, 0.01)
      dat.trial$ipsw.w.sv <- 1/s.score.sv
      tau.aipsw <- with(dat.trial, sum(ipsw.w.sv * ((A * (Y - Y1.t.sv) / sum(A * ipsw.w.sv) - ((1 - A) * (Y - Y0.t.sv)) / sum((1 - A) * ipsw.w.sv))))) +
        with(dat.rwe, mean(Y1.t.sv - Y0.t.sv))
    }
    ate$AIPSW <- tau.aipsw
  }

  if ("CW" %in% estimators) {
    tau.cw <- with(dat.trial, sum(q * ((A * (Y) / sum(A * q) - ((1 - A) * (Y)) / sum((1 - A) * q)))))
    ate$CW <- tau.cw
  }

  if ("ACW-t" %in% estimators) {
    if (sieve == FALSE) {
      tau.acw.t <- with(dat.trial, sum(q * ((A * (Y - Y1.t) / p1 - ((1 - A) * (Y - Y0.t)) / p0)))) + with(dat.rwe, mean(Y1.t - Y0.t))
    } else {
      tau.acw.t <- with(dat.trial, sum(q.t.sv * ((A * (Y - Y1.t.sv) / p1 - ((1 - A) * (Y - Y0.t.sv)) / p0)))) + with(dat.rwe, mean(Y1.t.sv - Y0.t.sv))
    }
    ate$ACW.t <- tau.acw.t
  }

  if ("ACW-b" %in% estimators) {
    if (sieve == FALSE) {
      tau.acw.b <- with(dat.trial, sum(q * ((A * (Y - Y1.b) / p1 - ((1 - A) * (Y - Y0.b)) / p0)))) + with(dat.rwe, mean(Y1.b - Y0.b))
    } else {
      tau.acw.b <- with(dat.trial, sum(q.b.sv * ((A * (Y - Y1.b.sv) / p1 - ((1 - A) * (Y - Y0.b.sv)) / p0)))) + with(dat.rwe, mean(Y1.b.sv - Y0.b.sv))
    }
    ate$ACW.b <- tau.acw.b
  }

  result <- c()
  result$ate <- ate
  result$hypers <- vector("list", 5)
  names(result$hypers) <- c("osel1.t", "osel0.t", "osel1.b", "osel0.b", "osel.ipsw")
  result$hypers$osel1.t <- osel1.t
  result$hypers$osel0.t <- osel0.t
  result$hypers$osel1.b <- osel1.b
  result$hypers$osel0.b <- osel0.b
  result$hypers$osel.ipsw <- osel.ipsw

  return(result)
}


