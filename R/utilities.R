expit <- function(x) {return (exp(x) / (1 + exp(x)))}

logit <- function(x) {return (log(x / (1 - x)))}

lamFun <- function(lam, moments, moments.bar) { ## vector lam and x
  qi <- (exp(moments %*% lam) / sum(exp(moments %*% lam)))[, 1]
  colSums(qi  * moments) - moments.bar
}

backward.sv <- function(osel1, osel0, beta1.sv, beta0.sv, gvec.trial, gvec.rwe){

  osel.all <- union(osel1[osel1 != 1], osel0[osel0 != 1])

  moms.t <- gvec.trial[, osel.all, drop = FALSE]
  moms <- gvec.rwe[, osel.all, drop = FALSE]
  moms.bar <- colMeans(moms)
  lam.hat <- searchZeros(matrix(rnorm(length(moms.bar) * 20, 0, 0.5), nrow = 20), lamFun, moments = moms.t, moments.bar = moms.bar)$x[1, ]
  if (!is.null(lam.hat)) {
    q.score <- exp(moms.t %*% lam.hat) / sum(exp(moms.t %*% lam.hat))
  } else {
    gvec.trial.new <- gvec.trial[, -1]
    gvec.rwe.new <- gvec.rwe[, -1]
    p.sv <- ncol(gvec.trial.new)

    beta0.new <- beta1.new <- numeric(p.sv + 1)
    beta1.new[osel1] <- abs(beta1.sv)
    beta0.new[osel0] <- abs(beta0.sv)
    beta.new <- pmax(beta1.new, beta0.new)[-1]

    while (length(osel.all) > 1 & is.null(lam.hat)) {
      min.beta <- min(beta.new[beta.new != 0])
      beta.new[beta.new == min.beta] <- 0
      osel.all <- which(beta.new != 0)
      moms.t <- gvec.trial.new[, osel.all, drop = FALSE]
      moms <- gvec.rwe.new[, osel.all, drop = FALSE]
      moms.bar <- colMeans(moms)
      lam.hat <- searchZeros(matrix(rnorm(length(moms.bar) * 20, 0, 0.5), nrow = 20), lamFun, moments = moms.t, moments.bar = moms.bar)$x[1, ]
    }

    if (!is.null(lam.hat)) {
      q.score <- exp(moms.t %*% lam.hat) / sum(exp(moms.t %*% lam.hat))
    } else {
      q.score <- 1 / nrow(gvec.trial.new)
    }
  }
  return(q.score)
}
