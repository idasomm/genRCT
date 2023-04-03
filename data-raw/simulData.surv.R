## code to prepare `simulData.surv` dataset goes here

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



'sample_prop' <- function(X, beta = 4.3) {
  expit(- beta - 0.5 * X[,1] - 0.5 * X[, 2] - 0.3 * X[, 3])
}


'trt_prop_RCT' <- function(X) {
  rep(0.5, nrow(X))
}

p <- 3

set.seed(12345)
X <- matrix(rnorm(300000 * p, 0, 1), ncol = p)
outlier <- data.frame(X) %>% mutate(id = c(1:300000)) %>%
  filter(X1 <= -1.645 | X1 >= 1.645) %>% filter(X2 <= -1.645 | X2 >= 1.645) %>%
  filter(X3 <= -1.645 | X3 >= 1.645) %>% dplyr::select(id)
X <- X[- outlier$id, ]
X <- data.frame(X) %>% filter(-4 <= X1 & X1 <= 4) %>% filter(-4 <= X2 & X2 <= 4) %>%
  filter(-4 <= X3 & X3 <= 4)
X <- X[1:2e5, ]
X <- as.matrix(X)


##### True survival curve
lambda.0 <- exp(- 3)
lambda.1 <- exp(- 3.7)

beta.1 <- c(-1, - 1, - 1.5)
beta.0 <- c(-1.8, - 1.5, - 1)

t.all <- seq(0.01, 50, by = 0.01)


##--------------------------------------------------------------------
## RCT data
X.pop <- X[1:50000, ]
eS <- sample_prop(X.pop, 4.3)
eS[eS >= 1] = 0.99
S <- sapply(eS, rbinom, n = 1, size = 1)
S.ind <- which(S == 1)
n.t <- length(S.ind)
X.t <- X[S.ind, ]
eA <- trt_prop_RCT(X.t) # treatment propensity score
A.t <- sapply(eA, rbinom, n = 1, size = 1)

## OS sample
n.p <- 5000
P.ind <- sample(50001:200000, size = n.p) ## RWD id
X.p <- X[P.ind, ]

Y.1 <- generate.survWeib.dep(rho = 1, X = X.t, lambda = lambda.1,  beta = beta.1,
                             XC = X.t, lambdaC = exp(- 4.5), betaC = c(-0.5, -1, -1), tmax = 400)
Y.0 <- generate.survWeib.dep(rho = 1, X = X.t, lambda = lambda.0,  beta = beta.0,
                             XC = X.t, lambdaC = exp(- 3.5), betaC = c(-0.5, -1, -1), tmax = 400)
Y.t <- Y.1$time * A.t + Y.0$time* (1 - A.t)
d.t <- Y.1$status * A.t + Y.0$status* (1 - A.t)


data.trial <- data.frame(Y = Y.t, d = d.t, A = A.t, X.t)
data.rwe <- data.frame(X.p)

simulData.surv <- vector("list", length = 2L)
names(simulData.surv) <- c("trial", "rwe")
simulData.surv$trial <- data.trial
simulData.surv$rwe <- data.rwe


usethis::use_data(simulData.surv, overwrite = TRUE)
