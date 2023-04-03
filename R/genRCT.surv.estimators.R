genRCT.surv.estimators <- function(Y.trial, d.trial, A.trial, X.trial, X.rwe, O.sel1 = NULL, O.sel0 = NULL,
                                     C.sel1 = NULL, C.sel0 = NULL, A.sel = NULL, eta = NULL, eta.vec = eta.vec, seed) {

  if(!is.null(seed)) set.seed(seed)
  n <- nrow(X.trial)
  p.trial <- ncol(X.trial)
  m <- nrow(X.rwe)
  n1 <- sum(A.trial == 1)
  n0 <- n - n1
  p1 <- n1 / n
  p0 <- 1 - p1

  dat.trial <- data.table(d = d.trial, A = A.trial, X = X.trial)
  dat.trial$Y <- Y.trial

  # Scale X
  X.all <- rbind(X.trial, X.rwe)
  X.all <- apply(X.all, 2, scale)
  X.trial <- X.all[1:n, ]
  X.rwe <- X.all[-c(1:n), ]

  # Sieve design matrix
  gvec.trial <- model.matrix(~ poly(X.trial, degree = 2, raw = TRUE))
  gvec.trial <- gvec.trial[, which(!duplicated(t(gvec.trial))), drop = FALSE]
  designX <- gvec.trial[, 2:ncol(gvec.trial)]
  gvec.rwe <- model.matrix(~ poly(X.rwe, degree = 2, raw = TRUE))
  gvec.rwe <- gvec.rwe[, which(!duplicated(t(gvec.rwe))), drop = FALSE]
  design.rwe <- gvec.rwe[, 2:ncol(gvec.rwe)]


  # Propensity score
  fit.ps <- glm(A.trial ~ X.trial, family = "binomial")
  ps <- as.vector(fit.ps$fitted)
  ps <- pmin(ps, 0.99)
  ps <- pmax(ps, 0.01)
  dat.trial$ps <- A.trial * ps + (1 - A.trial) * (1 - ps)

  # Sieve treatment propensity score
  if (is.null(A.sel)) {
    for(ns in 1:5) {
      fit.ps.sv <- cv.ncvreg(designX, y = A.trial, penalty = "SCAD", family = "binomial", nfolds = 5)
      beta.ps.sv <- as.numeric(coef(fit.ps.sv, s = 'lambda.min'))
      A.sel <- c(A.sel, union(1,which(beta.ps.sv != 0)))
    }
    tA <- table(A.sel)
    A.sel <- unique(sort(A.sel))[which(tA >= 3)]
  }

  fit.ps.sv <- glm(A.trial ~ cbind(1, designX)[, A.sel, drop = F] - 1, family = "binomial")
  ps.sv <- as.vector(fit.ps.sv$fitted)
  ps.sv <- pmin(ps.sv, 0.99)
  ps.sv <- pmax(ps.sv, 0.01)
  dat.trial$ps.sv <- A.trial * ps.sv + (1 - A.trial) * (1 - ps.sv)

  # IPSW
  X.all <- rbind(X.trial, X.rwe)
  delta <- c(rep(1, n), rep(0, m))
  #fit.ss <- glm(delta ~ X.all, family = 'binomial'(link = 'logit'))
  s.score <- as.vector(glm(delta ~ X.all, family = 'binomial'(link = 'logit'))$fitted)
  s.score <- pmin(s.score, 0.99)
  s.score <- pmax(s.score, 0.01)
  ipsw <- (1 - s.score[delta == 1]) / s.score[delta == 1]
  dat.trial$ipsw <- ipsw / sum(ipsw)

  # CW
  moment.bar <- colMeans(X.rwe)
  moms.t <- cbind(X.trial)# / (p1 * A.trial * ps + p0 * (1 - A.trial) * (1 - ps)))

  cnt <- 0
  while (cnt <= 50) {
    cnt <- cnt + 1
    lam.hat <- searchZeros(matrix(rnorm(length(moment.bar) * 20, 0, 0.25), nrow = 20), lamFun, moments = moms.t, moments.bar = moment.bar)$x[1,]
    if (!is.null(lam.hat)) break
  }

  if (is.null(lam.hat)) {
    #if (flag.boot == FALSE) warning('No lam.hat solutions')
    lam.hat <- rep(NA, p.trial)
  }
  q.score <- exp(X.trial %*% lam.hat) / sum(exp(X.trial %*% lam.hat))
  dat.trial$q <- q.score

  # Sieve q.score
  sv.X <- sieve.samp.2way(designX = designX, design.rwe = design.rwe, eta = eta, eta.vec = eta.vec)
  dat.trial$q.sv <- sv.X$q.sv
  eta <- sv.X$eta

  # Fit DR
  dr.1 <- est.dr(y = Y.trial[A.trial == 1], d = d.trial[A.trial == 1], x = X.trial[A.trial == 1, ], xnew = X.trial, ps = dat.trial$ps[A.trial == 1])
  dr.0 <- est.dr(y = Y.trial[A.trial == 0], d = d.trial[A.trial == 0], x = X.trial[A.trial == 0, ], xnew = X.trial, ps = dat.trial$ps[A.trial == 0])

  # Fit DACW
  dacw.1 <- est.DACW(y = Y.trial[A.trial == 1], d = d.trial[A.trial == 1], x = X.trial[A.trial == 1, ], xnew = X.rwe,
                     wt = dat.trial$q[A.trial == 1], ipsw = dat.trial$ipsw[A.trial == 1], ps = dat.trial$ps[A.trial == 1])
  dacw.0 <- est.DACW(y = Y.trial[A.trial == 0], d = d.trial[A.trial == 0], x = X.trial[A.trial == 0, ], xnew = X.rwe,
                     wt = dat.trial$q[A.trial == 0], ipsw = dat.trial$ipsw[A.trial == 0], ps = dat.trial$ps[A.trial == 0])


  # Fit DACW sieve
  dacw.sv.1 <- est.DACW.sv(y = Y.trial[A.trial == 1], d = d.trial[A.trial == 1], x = designX[A.trial == 1, ],
                           xnew = design.rwe, wt = dat.trial$q.sv[A.trial == 1], ps = dat.trial$ps.sv[A.trial == 1],
                           O.sel = O.sel1, C.sel = C.sel1)
  dacw.sv.0 <- est.DACW.sv(y = Y.trial[A.trial == 0], d = d.trial[A.trial == 0], x = designX[A.trial == 0, ],
                           xnew = design.rwe, wt = dat.trial$q.sv[A.trial == 0], ps = dat.trial$ps.sv[A.trial == 0],
                           O.sel = O.sel0, C.sel = C.sel0)

  return(list(surv1 = data.frame(time = dr.1$time, dr.1$surv, dacw.1$surv, dacw.sv.1$surv),
              surv0 = data.frame(time = dr.0$time, dr.0$surv, dacw.0$surv, dacw.sv.0$surv),
              hyper = list(O.sel1 = dacw.sv.1$O.sel, O.sel0 = dacw.sv.0$O.sel,
                           C.sel1 = dacw.sv.1$C.sel, C.sel0 = dacw.sv.0$C.sel, A.sel = A.sel, eta = eta)))
}
