
'est.dr' <- function(y, d, x, xnew, ps) {

  xnew <- as.matrix(xnew)

  # Fit coxph
  fit <- coxph(Surv(y, d) ~ x)
  beta <- fit$coefficients
  risk <- c(exp(x %*% beta))
  risk.new <- c(exp(xnew %*% beta))

  base <- basehaz(fit, centered = FALSE)
  time <- base$time
  h0 <- base$hazard
  dh0 <- c(h0[1], diff(h0))
  dhazard <- outer(risk, dh0, "*")
  dhazard.new <- outer(risk.new, dh0, "*")
  hazard <- outer(risk, h0, "*")
  hazard.new <- outer(risk.new, h0, "*")

  # Fit probability of censoring
  fitC <- coxph(Surv(y, 1 - d) ~ x)
  betaC <- fitC$coefficients
  riskC <- c(exp(x %*% betaC))

  baseC <- basehaz(fitC, centered = FALSE)
  h0C <- baseC$hazard
  dh0C <- c(h0C[1], diff(h0C))
  dhazardC <- outer(riskC, dh0C, "*")
  hazardC <- outer(riskC, h0C, "*")

  ## Martingale
  # Counting process for survival time
  y.ai <- outer(y, time, ">=")
  dy.ai <- outer(y, time, "==")
  dn.ai  <- dy.ai * d

  # Counting process for censoring time
  dyC.ai <- outer(y, baseC$time, "==")
  dnC.ai  <- dyC.ai * (1 - d)
  dmC.ai <- dnC.ai - dhazardC * y.ai

  # Denominator of the robust estimator
  denom1 <- exp(- hazard.new)
  denom2 <- 1 / ps * y.ai * exp(hazardC)
  denom3 <- 1 / ps * exp(- hazard)
  argC <- t(apply(exp(hazard + hazardC) * dmC.ai, 1, cumsum))
  denom4 <- denom3 * argC

  n <- length(xnew[, 1])
  denom <- colMeans(denom1) + colSums(denom2 - denom3 + denom4, na.rm = TRUE) / n
  denom <- pmax(pmin(denom, 0.99), 0.01)

  # Numerator of the robust estimator
  num1 <- denom1 * dhazard.new
  num2 <- 1 / ps * dn.ai * exp(hazardC)
  num3 <- denom3 * dhazard
  num4 <- num3 * argC
  num <- colMeans(num1) + colSums(num2 - num3 + num4, na.rm = TRUE) / n

  hazard.dr <- cumsum(num / denom)
  surv.dr <- exp(- hazard.dr)

  # Naive estimator
  denom.nv <- colSums(denom2, na.rm = T) / n
  num.nv <- colSums(num2, na.rm = T) / n
  surv.nv <- exp(- cumsum(num.nv / denom.nv))

  #time.evt <- coxph.detail(fit)$time
  #t.idx <- time %in% time.evt
  out <- data.frame(surv.nv = surv.nv, denom.nv = denom.nv, surv.dr = surv.dr, denom.dr = denom)
  #out[is.nan(as.matrix(out))] <- NA
  #out <- fill(out)
  out <- apply(out, 2, function(x) pmax(x, 0))
  out <- apply(out, 2, function(x) pmin(x, 1))
  out <- apply(out, 2, function(x) mono.dec(x))

  return(list(time = time, surv = out))
}

'est.DACW' <- function(y, d, x, xnew, wt, ipsw, ps) {

  # Fit coxph
  fit <- coxph(Surv(y, d) ~ x)
  beta <- fit$coefficients
  risk <- c(exp(x %*% beta))
  risk.new <- c(exp(xnew %*% beta))

  base <- basehaz(fit, centered = FALSE)
  time <- base$time
  h0 <- base$hazard
  dh0 <- c(h0[1], diff(h0))
  dhazard <- outer(risk, dh0, "*")
  dhazard.new <- outer(risk.new, dh0, "*")
  hazard <- outer(risk, h0, "*")
  hazard.new <- outer(risk.new, h0, "*")

  # Fit probability of censoring
  fitC <- coxph(Surv(y, 1 - d) ~ x)
  betaC <- fitC$coefficients
  riskC <- c(exp(x %*% betaC))

  baseC <- basehaz(fitC, centered = FALSE)
  h0C <- baseC$hazard
  dh0C <- c(h0C[1], diff(h0C))
  dhazardC <- outer(riskC, dh0C, "*")
  hazardC <- outer(riskC, h0C, "*")

  index_to_grid <- unlist(lapply(y, function(x) sum(x >= base$time)))
  wC <- exp(h0C[index_to_grid] * riskC)

  ## Martingale
  # Counting process for survival time
  y.ai <- outer(y, time, ">=")
  dy.ai <- outer(y, time, "==")
  dn.ai  <- dy.ai * d

  # Counting process for censoring time
  dyC.ai <- outer(y, baseC$time, "==")
  dnC.ai  <- dyC.ai * (1 - d)
  dmC.ai <- dnC.ai - dhazardC * y.ai

  # Normalize weights
  wt.norm <- wt / ps
  wt.norm <- wt.norm / sum(wt.norm)

  ipsw.norm <- ipsw / ps
  ipsw.norm <- ipsw.norm / sum(ipsw.norm)

  # Denominator of the robust estimator
  denom1 <- wt.norm * exp(hazardC) * y.ai
  denom2 <- wt.norm * exp(- hazard)
  denom3 <- exp(- hazard.new)
  argC <- t(apply(exp(hazard + hazardC) * dmC.ai, 1, cumsum))
  denom4 <- denom2 * argC

  # plot(t.all, surv.true$surv1, type='l', ylim=c(0,1))
  # lines(time, colSums(denom1), col = 2)
  # lines(time, colSums(denom2), col = 3)
  # lines(time, colMeans(denom3), col = 4)
  # lines(time, colSums(denom4), col = 5)
  #

  num1 <- wt.norm * exp(hazardC) * dn.ai
  num2 <- denom2 * dhazard
  num3 <- denom3 * dhazard.new
  num4 <- num2 * argC

  denom <- colSums(denom1 - denom2 + denom4, na.rm = TRUE) + colMeans(denom3, na.rm = T)
  denom <- pmax(pmin(denom, 0.99), 0.01)
  num <- colSums(num1 - num2 + num4, na.rm = TRUE) + colMeans(num3, na.rm = T)
  hazard.dacw <- cumsum(num / denom)
  surv.dacw <- exp(- hazard.dacw)

  # IPSW estimator
  denom.ipsw <- colSums(ipsw.norm * exp(hazardC) * y.ai, na.rm = TRUE)
  num.ipsw <-  colSums(ipsw.norm * exp(hazardC) * dn.ai, na.rm = TRUE)
  surv.ipsw <- exp(- cumsum(num.ipsw / denom.ipsw))

  # CW estimator
  denom.cw <- colSums(denom1, na.rm = TRUE)
  num.cw <- colSums(num1, na.rm = TRUE)
  surv.cw <- exp( - cumsum(num.cw / denom.cw))

  # OR estimator
  denom.or <- colMeans(denom3, na.rm = TRUE)
  num.or <- colMeans(num3, na.rm = TRUE)
  surv.or <- exp( - cumsum(num.or / denom.or))

  #time.evt <- coxph.detail(fit)$time
  #t.idx <- base$time %in% time.evt
  out <- data.frame(surv.ipsw = surv.ipsw, denom.ipsw = denom.ipsw,
                    surv.cw = surv.cw, denom.cw = denom.cw,
                    surv.or = surv.or, denom.or = denom.or,
                    surv.dacw = surv.dacw, denom.dacw = denom)
  #out[is.nan(as.matrix(out))] <- NA
  #out <- fill(out)
  out <- apply(out, 2, function(x) pmax(x, 0))
  out <- apply(out, 2, function(x) pmin(x, 1))
  out <- apply(out, 2, function(x) mono.dec(x))

  return(list(time = time, surv = out))
}


'est.DACW.sv' <- function(y, d, x, xnew, wt, ps, O.sel = NULL, C.sel = NULL) {

  if (is.null(O.sel) & is.null(C.sel)) {
    for(ns in 1:5) {
      fit <- cv.ncvsurv(x, y = Surv(y, d), penalty = "SCAD", nfolds = 5)
      beta <- as.numeric(coef(fit, s = 'lambda.min'))
      O.sel <- union(O.sel, which(beta != 0))
      fitC <- cv.ncvsurv(x, y = Surv(y, 1 - d), penalty = "SCAD", nfolds = 5)
      betaC <- as.numeric(coef(fitC, s = 'lambda.min'))
      C.sel <- union(C.sel, which(betaC != 0))
    }
    #tO <- table(O.sel)
    #tC <- table(C.sel)
    #O.sel <- unique(sort(O.sel))[which(tO >= 3)]
    #C.sel <- unique(sort(C.sel))[which(tC >= 3)]
  }



  #Refit the selected basis for the Outcome coxph

  if (length(O.sel) == 0) {
    fit <- coxph(Surv(y, d) ~ 1)
    risk <- rep(1, nrow(x))
    risk.new <- rep(1, nrow(xnew))
  } else {
    fit <- coxph(Surv(y, d) ~ x[, O.sel, drop = F])
    beta <- fit$coefficients
    risk <- c(exp(x[, O.sel, drop = F] %*% beta))
    risk.new <- c(exp(xnew[, O.sel, drop = F] %*% beta))
  }


  base <- basehaz(fit, centered = FALSE)
  time <- base$time
  h0 <- base$hazard
  dh0 <- c(h0[1], diff(h0))
  dhazard <- outer(risk, dh0, "*")
  dhazard.new <- outer(risk.new, dh0, "*")
  hazard <- outer(risk, h0, "*")
  hazard.new <- outer(risk.new, h0, "*")

  # ReFit probability of censoring
  if (length(C.sel) == 0) {
    fitC <- coxph(Surv(y, 1 - d) ~ 1)
    riskC <- rep(1, nrow(x))
  } else {
    fitC <- coxph(Surv(y, 1 - d) ~ x[, C.sel, drop = F])
    betaC <- fitC$coefficients
    riskC <- c(exp(x[, C.sel, drop = F] %*% betaC))
  }


  baseC <- basehaz(fitC, centered = FALSE)
  h0C <- baseC$hazard
  dh0C <- c(h0C[1], diff(h0C))
  dhazardC <- outer(riskC, dh0C, "*")
  hazardC <- outer(riskC, h0C, "*")

  index_to_grid <- unlist(lapply(y, function(x) sum(x >= base$time)))
  wC <- exp(h0C[index_to_grid] * riskC)

  ## Martingale
  # Counting process for survival time
  y.ai <- outer(y, time, ">=")
  dy.ai <- outer(y, time, "==")
  dn.ai  <- dy.ai * d

  # Counting process for censoring time
  dyC.ai <- outer(y, baseC$time, "==")
  dnC.ai  <- dyC.ai * (1 - d)
  dmC.ai <- dnC.ai - dhazardC * y.ai

  # Normalize weights
  wt.norm <- wt / ps
  wt.norm <- wt.norm / sum(wt.norm)


  # Denominator of the robust estimator
  denom1 <- wt.norm * exp(hazardC) * y.ai
  denom2 <- wt.norm * exp(- hazard)
  denom3 <- exp(- hazard.new)
  argC <- t(apply(exp(hazard + hazardC) * dmC.ai, 1, cumsum))
  denom4 <- denom2 * argC
  #

  num1 <- wt.norm * exp(hazardC) * dn.ai
  num2 <- denom2 * dhazard
  num3 <- denom3 * dhazard.new
  num4 <- num2 * argC

  #denom <- mono.dec(colSums(denom1, na.rm = TRUE)) - colSums(denom2 - denom4, na.rm = TRUE) + mono.dec(colMeans(denom3, na.rm = TRUE))
  denom <- colSums(denom1 - denom2 + denom4, na.rm = TRUE) + colMeans(denom3, na.rm = T)
  #denom <- pmax(pmin(denom, 0.99), 0.01)
  #denom <- mono.dec(denom)
  num <- colSums(num1 - num2 + num4, na.rm = TRUE) + colMeans(num3, na.rm = T)
  hazard.dacw.sv <- cumsum(num / denom)
  surv.dacw.sv <- mono.dec(exp(- hazard.dacw.sv))

  #time.evt <- coxph.detail(fit)$time
  #t.idx <- base$time %in% time.evt
  out <- data.frame(surv.dacw.sv = surv.dacw.sv, denom.dacw.sv = denom)
  #out[is.nan(as.matrix(out))] <- NA
  #out <- fill(out)
  out <- apply(out, 2, function(x) pmax(x, 0.001))
  out <- apply(out, 2, function(x) pmin(x, 1))

  # plot(t.all, surv.true$surv1, type='n', ylim=c(0,1))
  # lines(time[t.idx], mono.dec(colSums(denom1, na.rm=T))[t.idx], col = 2)
  # lines(time[t.idx], mono.dec(colSums(denom2, na.rm=T))[t.idx], col = 3)
  # lines(time[t.idx], mono.dec(colMeans(denom3,na.rm=T))[t.idx], col = 4)
  # lines(time, colSums(denom4), col = 5)

  out <- apply(out, 2, function(x) mono.dec(x))

  return(list(time = time, surv = out, O.sel = O.sel, C.sel = C.sel))
  #return(list(time = time, surv = out, O.sel = O.sel, C.sel = C.sel))

}




'generate.survWeib.dep' <- function(rho, X, lambda, beta, XC, lambdaC, betaC, tmax) {

  if(!is.matrix(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  u <- runif(n)
  Tlat <- (-log(u) / (lambda * exp(X %*% beta)))^(1/rho)

  # censoring times
  #C <- rexp(n) * exp(- rateC)
  #C <- rexp(n, rate = rateC)
  u <- runif(n)
  C <- (-log(u) / (lambdaC * exp(XC %*% betaC)))^(1/rho)
  C <- pmin(C, tmax)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  cat('# Proportion of censoring','\n')
  print(mean(1 - status))


  # data set
  return(data.frame(time = time, status = status, Tlat = Tlat, C = C))
}


'generate.survAFT' <- function(X, beta, rateC, tmax) {

  if(!is.matrix(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  Tlat <- exp(X %*% beta + rnorm(n))

  # censoring times
  #C <- rexp(n) * exp(- rateC)
  #C <- rexp(n, rate = rateC)
  u <- runif(n)
  C <- -log(u) / exp(-rateC)
  C <- pmin(C, tmax)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  cat('# Proportion of censoring','\n')
  print(mean(1 - status))

  # data set
  return(list(time = time, status = status))
}


'get.rmst' <- function(time, surv, tau){

  time <- round(time, 5)
  tau <- round(tau, 5)
  L <- length(tau)
  rmst <- rep(NA, L)
  data <- data.frame(time = time, surv = surv) %>% arrange(time)
  for(l in 1:L) {
    trunc.data <- data %>% filter(time <= tau[l])
    rmst[l] <- sum(diff(c(0, trunc.data$time, tau[l])) * c(1, trunc.data$surv))
  }

  sum(diff(c(0, trunc.data$time, tau[l])) * c(1, trunc.data$surv))
  sum(diff(c(0, trunc.data$time)) * c(trunc.data$surv))


  names(rmst) <- paste0("tau = ", tau)
  return(rmst)
}

'get.rmst.est' <- function(fit.obj, est, tau) {
  rmst1 <- get.rmst(time = fit.obj$surv1$time, surv = fit.obj$surv1[[est]], tau = tau)
  rmst0 <- get.rmst(time = fit.obj$surv0$time, surv = fit.obj$surv0[[est]], tau = tau)
  rmst <- rbind(rmst1, rmst0, rmst1 - rmst0)
  rownames(rmst) <- c("rmst.1", "rmst.0", "rmst.diff")
  return(rmst)
}


'rmst_to_save' <- function(fit.obj, tau) {

  rmst.est <- vector("list", length = 14)
  est.names <- names(fit.obj$surv1)[- 1]
  names(rmst.est) <- est.names

  for(name in est.names) {
    rmst.est[[name]] <- as.vector(get.rmst.est(fit.obj = fit.obj, est = name, tau = tau))
  }

  rmst <- matrix(unlist(rmst.est), ncol = 14)

  colnames(rmst) <- est.names
  name1 <- rep(paste0('t', tau), each = 3)
  name2 <- rep(c("1","0","diff"), 5)
  rownames(rmst) <- paste(name1, name2, sep='.')
  return(rmst)
}


'mono.dec' <- function(x) {
  x <- as.vector(x)
  xnew <- c(1, x)
  lower.mat <- outer(xnew, xnew, FUN="pmin")
  lower.mat[upper.tri(lower.mat, diag = F)] <- NA
  out <- apply(lower.mat, 1, min, na.rm = T)
  return(out[-1])
}


'sieve.samp.2way' <- function(designX, design.rwe, par0 = NULL, eta = NULL, eta.vec = c(0,0.001,0.005)){

  moms.bar <- colMeans(design.rwe)
  moms.t <- designX

  sieve <- sieve.samp.q2(moms.t = moms.t, moms.bar = moms.bar, eta = eta, eta.vec = eta.vec)


  return(list(q.sv = sieve$q2, eta = sieve$eta))
}


'sieve.samp.q2' <- function(moms.t, moms.bar, eta = NULL, eta.vec){

  par0 <- nleqslv::searchZeros(matrix(rnorm(length(moms.bar)*10,0,0.1),nrow = 10),lamFun,
                               moments = moms.t, moments.bar = moms.bar)$x[1,]

  if(is.null(eta)){
    #cvObj <- cvPEE(par0, moms.t = moms.t, moms.bar = moms.bar, eta.vec = eta.vec, nfolds = 5)
    cvObj <- cvPEE2(par0, moms.t = moms.t ,moms.bar = moms.bar, eta.vec = eta.vec, train.prop = 1/2, cv.n = 10)
    eta <- cvObj$eta.min
  }

  lam.hat0 <- matrix(NA, nrow = 5, ncol = ncol(moms.t))
  par0 <- matrix(rnorm(ncol(moms.t)*5,0,0.1), nrow = 5)

  #S.sel <- NULL
  for(i in 1:nrow(par0)){
    lam.hat0[i,] <- solveEE(par0[i,], moms.t = moms.t,moms.bar = moms.bar,eta1 = eta)
    #S.sel <- c(S.sel, which(lam.hat0[i,] != 0))
  }

  # tS <- table(S.sel)
  # S.sel <- unique(sort(S.sel))[which(tS >= 3)]
  #

  lam.hat0 <- round(lam.hat0, 5L)
  notdups <- !duplicated(lam.hat0)
  lam.hat1 <- lam.hat0[notdups, , drop = F]

  if (nrow(lam.hat1) > 1) {
    dd <- as.matrix(dist(lam.hat0, diag = T, upper = T))
    d.min <- which.min(colSums(dd)[notdups])
    lam.hat1 <- lam.hat1[d.min, ]
  }

  lam.hat1 <- c(lam.hat1)
  S.sel <- which(lam.hat1!=0)

  if(length(S.sel) <= 2){
    sel <- 1:length(par0)
  }

  if(length(S.sel)!=ncol(moms.t)){
    lam.hat00 <- lam.hat1[S.sel]
    moms.t1 <- moms.t[,S.sel]
    moms.bar1 <- moms.bar[S.sel]

    repeat{

      ## Re-calibarate
      lam.hat <- nleqslv::searchZeros(matrix(rnorm(length(lam.hat00)*10,0,0.1),nrow = 10),lamFun,
                                      moments = moms.t1, moments.bar = moms.bar1)$x[1,]
      if(!is.null(lam.hat)) {q2 <- exp(moms.t1%*%lam.hat)/sum(exp(moms.t1%*%lam.hat)); break}

      # Delete the smallest lambda among the selected
      lam.min <- which.min(abs(lam.hat00))
      lam.hat00 <-lam.hat00[-lam.min]
      moms.t1 <- moms.t1[,-lam.min]
      moms.bar1 <- moms.bar1[-lam.min]

      if(length(lam.hat00) == 2) {
        q2 <- exp(moms.t%*%lam.hat1)/sum(exp(moms.t%*%lam.hat1)); break
      }
    }

  } else {q2 <- exp(moms.t%*%lam.hat1)/sum(exp(moms.t%*%lam.hat1))}

  return(list(q2 = q2, eta = eta))
}


solveEE <- function(par0,moms.t,moms.bar,eta1){
  nn <- nrow(moms.t)
  pp <- ncol(moms.t)
  En <- matrix(0,pp,pp)
  cnt1 <- 1
  cnt2 <- 0
  epsilon <- 10^(-6)
  par <- par0
  par.var <- 0.1

  while(cnt1 <= 100){
    cnt1 <- cnt1 + 1
    cnt2 <- cnt2 + 1

    if(cnt2 == 200){
      par.var <- par.var + 0.05
      cnt2 <- 0
    }
    if(par.var > 1){cat('SolveEE convergence not met','\n'); break}
    ee0 <- UEE(par0,moments = moms.t,moments.bar = moms.bar)
    jac <- -lamJac(par0,moments = moms.t,moments.bar = moms.bar)
    diag(En) <- scadPen(par0,eta = eta1, a = 3.7)/(epsilon+abs(par0))
    #diag(En) <- eta1/(epsilon+abs(par0))

    #((abs(par0)<eta1) + (abs(lam)>=eta)*((a*eta) > abs(lam))*(a*eta - abs(lam))/(a-1))

    #diag(En) <- ridgePen(par0,eta = eta1)/(epsilon+abs(par0))
    res <- try(ginv(jac+En),silent = TRUE)

    #res <- ginv(jac+En)
    if(inherits(res, "try-error") )
    { #error handling code
      par0 <- rnorm(length(par0),0,par.var)
      cnt1 <- ceiling(cnt1/2)+1
      next;
      # break;
    }
    upd <- (res%*%(ee0-En%*%par0))[,1]
    # print(round(upd,3))
    par <- par0 + upd
    # if(any(abs(par) > 50)){print("bad");break;}
    if( sum( abs(par-par0)  )<10^(-5)){
      par <- par*(abs(par) >= 0.001)
      # print(norm(jac,'2'));
      break;
    }
    if( sum( abs(par-par0)  )>100 ) { #error handling code
      par0 <- rnorm(length(par0),0,par.var)
      cnt1 <- ceiling(cnt1/2)+1
      next;
    }
    par0 <- par

  }


  # print(cnt1)
  par[which( abs(par) <10^(-3))]<-0
  return(par)
}


tunePEE <- function(par0,moms.t,moms.bar,eta.vec,FLAG=0){
  nn <- nrow(moms.t)
  ll <- length(eta.vec)
  value.vec <- numeric(ll)
  lam.mat <- matrix(0,nrow = ll,ncol = length(par0))
  for(j in 1:length(eta.vec)){
    eta1 <- eta.vec[j]
    lam.mat[j,] <- solveEE(par0,moms.t = moms.t, moms.bar = moms.bar,eta1 = eta1)
    # check <- lamEE(lam = lam.mat[j,],moments = moms.t,moments.bar = moms.bar,eta = eta1)
    # if(sum(check^2)>100){
    #   # print("bad")
    #   lam.init <- par0 + rnorm((length(moms.bar)),0,0.1)
    #   lam.mat[j,] <- solveEE(lam.init,moms.t = moms.t, moms.bar = moms.bar,eta1 = eta1)
    #
    #   lam.hat1 <- searchZeros(lam.init, lamEE, moments = moms.t1, moments.bar = moms.bar1,eta=0,
    #                           control=list(trace=0))$x[1,]
    #   ct = 1
    #   while(is.null(lam.hat1) & ct <= 10){
    #     lam.init <- par0 + matrix(rnorm((length(moms.bar))*10,0,0.2),nrow = 10)
    #     lam.hat1 <- searchZeros(lam.init, lamEE, moments = moms.t1, moments.bar = moms.bar1,eta=0,
    #                             control=list(trace=0))$x[1,]
    #     ct <- ct +1
    #     print(ct)
    #   }
    #   if(!is.null(lam.hat1)){
    #     lam.mat[j,] <- lam.hat1
    #   }else{FLAG = 1}
    # }
    value.vec[j] <- lossFun(lam.mat[j,],moments = moms.t,moments.bar = moms.bar)
  }

  min.ind <- which.min(value.vec)
  eta.min <- eta.vec[min.ind]
  lam.hat <- lam.mat[min.ind,]
  return(list(lam.hat = lam.hat,eta.min = eta.min, loss.eta = value.vec, eta.vec = eta.vec))
}

UEE <- function(lam, moments, moments.bar){
  ## Estimating equation for lambda
  # qi <- (exp(moments%*%lam)/sum(exp(moments%*%lam)))[,1]
  # ee <- colSums(qi*moments) - moments.bar
  # ee
  q <- exp(moments%*%lam)[,1]
  tmp <- t(apply(moments,1,function(x){x - moments.bar}))
  ee <- colSums(diag(q)%*%tmp)/sum(q)
  ee
}

lamJac <-function(lam,moments,moments.bar){
  ## Derivative of lamEE
  q <- exp(moments%*%lam)[,1]
  tmp <- t(apply(moments,1,function(x){x - moments.bar}))
  pp <- ncol(moments)
  jac <- matrix(0,pp,pp)
  for(i in 1:nrow(moments)){
    jac <- jac + q[i]*tmp[i,]%*%t(moments[i,])
  }
  jac/sum(q)
}

scadPen <- function(lam, eta, a = 3.7){
  penaltyd <- ((abs(lam)<eta) + (abs(lam)>=eta)*((a*eta) > abs(lam))*(a*eta - abs(lam))/(a-1))
  eta*penaltyd
}

lossFun <- function(lam, moments, moments.bar){
  q <- exp(moments%*%lam)[,1]
  tmp <- cbind(q,moments)
  tmp1 <- apply(tmp,1,function(x){as.numeric(x[1]*crossprod((x[2:length(x)] - moments.bar)))})
  loss <- sum(tmp1)/sum(q)
  loss
}


Create_Folds <- function(y, k = 10, list = TRUE, returnTrain = FALSE)
{
  if (class(y)[1] == "Surv")
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}


balance.folds <- function(x, train.prop = 3/4, error.max = 0.1, max.iter = 5000) {

  sdx <- matrixStats::colSds(x)
  n <- nrow(x)
  n.train <- round(n * train.prop, digits = 0)
  # n1 <- sum(trt)
  # n0 <- sum(1 - trt)
  # m1 <- round(n1 * train.prop, digits = 0)
  # m0 <- round(n0 * train.prop, digits = 0)
  # n <- n1 + n0
  # m <- m1 + m0
  # id1 <- (1:n)[trt == 1]
  # id0 <- (1:n)[trt == 0]
  error <-  errorbest <- Inf
  bestid.test <- rep(NA, n - n.train)
  iter <- 0

  while ((error > error.max | is.na(error) == TRUE) & iter <= max.iter) {
    id.test <- c(sample(x = 1:n, size = n - n.train, replace = FALSE))
    x.test <- x[id.test, , drop = FALSE]
    x.train <- x[-id.test, , drop = FALSE]

    diffx <- colMeans(x.train) - colMeans(x.test)
    error <- max(abs(diffx / sdx))

    if (is.na(error) == FALSE & error < errorbest) {
      bestid.test <- id.test
      errorbest <- error
    }
    iter <- iter + 1
  }

  if (all(is.na(bestid.test)) == TRUE) stop("No balanced data split found.")

  if(iter == max.iter + 1){

    x.test <- x[bestid.test, , drop = FALSE]
    x.train <- x[-bestid.test, , drop = FALSE]

    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(x.train = x.train, x.test = x.test))
}

balance.boots <- function(x, A.trial, trt = 1, error.max = 0.1, max.iter = 5000, rwe = FALSE) {

  n <- nrow(x)
  if (rwe == FALSE) {
    xnew <- x[A.trial == trt,]
  } else {
    xnew <- x
  }
  sdx <- matrixStats::colSds(xnew)


  error <-  errorbest <- Inf
  bestid.samp <- rep(NA, n)
  iter <- 0

  while ((error > error.max | is.na(error) == TRUE) & iter <= max.iter) {
    if (rwe == FALSE) {
      id.samp <- sample(which(A.trial == trt), replace = TRUE)
    } else {
      id.samp <- sample(1:n, replace = TRUE)
    }
    x.samp <- x[id.samp, , drop = FALSE]

    diffx <- colMeans(x.samp) - colMeans(xnew)
    error <- max(abs(diffx / sdx))

    if (is.na(error) == FALSE & error < errorbest) {
      bestid.samp <- id.samp
      errorbest <- error
    }
    iter <- iter + 1
  }

  if (all(is.na(bestid.samp)) == TRUE) stop("No balanced data split found.")

  if(iter == max.iter + 1){

    id.samp <- bestid.samp

    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(id.samp)
}


cvPEE2 <- function(par0,moms.t,moms.bar,eta.vec, train.prop = 1/2, cv.n = 10){

  value_lam <- function(eta1){
    value.vec <- numeric(cv.n)
    #nk <- numeric(cv.n)
    for(cv.k in 1:cv.n){
      split.x <- balance.folds(moms.t, train.prop = train.prop)
      moms.train <- split.x$x.train
      moms.test <- split.x$x.test

      lam.hat <- solveEE(par0,moms.t = moms.train,moms.bar = moms.bar,eta1 = eta1)
      #nk[cv.k] <- length(lam.hat != 0)
      value.vec[cv.k] <- lossFun(lam.hat,moments = moms.test,moments.bar = moms.bar)
    }
    return(c(median(value.vec)))
  }
  cv.value <- sapply(eta.vec,value_lam)
  min.ind <- which.min(cv.value)
  eta.min <- eta.vec[min.ind]
  lam.hat <- solveEE(par0,moms.t = moms.t, moms.bar = moms.bar,eta1 = eta.min)
  return(list(lam.hat = lam.hat,eta.min = eta.min, loss.eta = cv.value, eta.vec = eta.vec))
}


cvPEE <- function(par0,moms.t,moms.bar,eta.vec,nfolds = 2){

  nn <- nrow(moms.t)
  CV.ind <- Create_Folds(1:nn,k=nfolds)

  value_lam <- function(eta1){
    value.vec <- numeric(nfolds)
    k = 1
    while(k <= nfolds){
      test.ind <- CV.ind[[k]]
      train.ind <- setdiff(1:nn,test.ind)
      moms.train <- moms.t[train.ind,]
      lam.hat <- solveEE(par0,moms.t = moms.train,moms.bar = moms.bar,eta1 = eta1)
      value.vec[k] <- lossFun(lam.hat,moments = moms.t[test.ind,],moments.bar = moms.bar)
      k <- k + 1

    }
    return(median(value.vec))
  }
  cv.value <- sapply(eta.vec,value_lam)
  min.ind <- which.min(cv.value)
  eta.min <- eta.vec[min.ind]
  lam.hat <- solveEE(par0,moms.t = moms.t, moms.bar = moms.bar,eta1 = eta.min)
  return(list(lam.hat = lam.hat,eta.min = eta.min, loss.eta = cv.value, eta.vec = eta.vec))
}

alassoPen <-function(weights,eta=0){
  ## ALASSO penalty derivative
  eta*abs(weights)
}

ridgePen <- function(lam, eta=0){
  penaltyd <- lam#/abs(lam_0)
  eta*penaltyd
}
