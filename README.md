
# genRCT

### Generalizing the average treatment effect (ATE) from the trial using observational studies (OS)

### Installation

``` r
devtools::install_github("idasomm/genRCT")
```

### Usage

genRCT(Y.trial, X.trial, A.trial, Y.rwe, X.rwe, A.rwe, family =
“gaussian”, estimators = c(“Naive”, “IPSW”, “AIPSW”, “CW”, “ACW-t”,
“ACW-b”), sieve = TRUE, inference = TRUE, n.boot = 100, conf.level =
0.05, seed = NULL, plot.boot = TRUE, verbose = TRUE)

### Help functions

``` r
library(genRCT)
?genRCT
```

### Arguments

| Argument   |                                                                                                                                                                                                                                        |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Y.trial    | Observed outcome from a trial; vector of size \(n\) (the trial sample size).                                                                                                                                                           |
| A.trial    | Treatment received from a trial; vector of size \(n\).                                                                                                                                                                                 |
| X.trial    | Matrix of \(p\) baseline covariates from a trial; dimension \(n\) by \(p\).                                                                                                                                                            |
| Y.rwe      | Observed outcome from OS; if obtained, vector of size \(m\) (OS sample size); otherwise, set Y.rwe = NULL.                                                                                                                             |
| A.rwe      | Treatment received from OS; if obtained, vector of size \(m\); otherwise, set Y.rwe = NULL.                                                                                                                                            |
| X.rwe      | Matrix of \(p\) baseline covariates from OS; dimension $m by \(p\).                                                                                                                                                                    |
| family     | The type of outcome; “gaussian” for continuous outcome or “binomial” for binary outcome. Default is “gaussian”.                                                                                                                        |
| estimators | A vector of one or multiple methods to estimate the ATE. Allowed values are “Naive”, “IPSW”, “AIPSW”, “CW”, “ACW-t”, “ACW-b”. The “ACW-b” is allowed only when both “Y.rwe” and “A.rwe” are obtained. Default specifies all 6 methods. |
| sieve      | A logical value indicating whether the method of sieves are used for estimating sampling score and outcome models. Used only if estimators = “AIPSW or”ACW-t" or “ACW-b”. Default is TRUE.                                             |
| inference  | A logical value indicating whether inference for the ATE via bootstrap should be provided. Default it TRUE.                                                                                                                            |
| n.boot     | A numeric value indicating the number of bootstrap samples used. This is only relevant if inference = TRUE. Default is 100.                                                                                                            |
| conf.level | The level of bootstrap confidence interval; Default is 0.05.                                                                                                                                                                           |
| seed       | An optional integer specifying an initial randomization seed for reproducibility. Default is NULL, corresponding to no seed.                                                                                                           |
| plot.boot  | A logical value indicating whether histograms of the bootstrap samples should be produced. Default is TRUE.                                                                                                                            |
| verbose    | A logical value indicating whether intermediate progress messages should be printed. Default is TRUE.                                                                                                                                  |

### Value

| Value |                                                                                                                                                                |
| ----- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| fit   | A table of estimated ATEs with bootstrap SE and confidence interval.                                                                                           |
| plot  | A set of histograms displaying the distribution of the bootstrapped estimates. The red vertical reference lines represent the estimated ATEs from each method. |

### Example

``` 

library(genRCT)
library(data.table)
library(nleqslv)
library(ncvreg)
library(MASS)
library(Matrix)
library(kableExtra)
library(tidyverse)

# Data generating functions

'outcome_mod_lognorm' <- function(X,A,sigma = 0.25){
  n <- length(A)
  return( -100 + 27.4*X[,3]*A + 13.7*X[,4] + 5*X[,4]*A + exp(sigma*rnorm(n)))
}

'outcome_mod' <- function(X,A){
  expit(1 - 2*A*X[,3] - X[,3] - 0.5*A*X[,4] + X[,4])
}

'sample_prop' <- function(X,beta = 2){
  #tmp <- -beta + 0.5*X[,1] + 0.3*X[,2] - 0.5*X[,3] -0.4* X[,4]
  exp(-beta + X[,1] + 0.6*X[,2] - 0.5*X[,3])
}


'trt_prop' <- function(X){
  expit(-X[,1] + 0.4*X[,2] - 0.25*X[,3] - 0.1*X[,4])
}


'gen_Z' <- function(X){
  z1 <- exp(X[,1]/2)
  z2 <- (X[,1]+X[,4]+20)^2
  z3 <- X[,2]/(1+2*exp(X[,3]))
  #z3 <- (X[,1]*X[,3]/2 + 2)^3
  z4 <- (X[,2]+X[,5]+20)^2
  return(cbind(z1,z2,z3,z4))
}



######################
#  gaussian outcome  #
######################


p <- 4
n.p <- 2000
n <- 2e4 ## population size
n2 <-2e5 ## population size

beta0 <-5.3

## Generating data
set.seed(12345)
X <- matrix(rnorm(n*p,1,1),ncol = p)
X <- scale(X) + 1
eS <- sample_prop(X, beta0)

#print(sum(eS>=1))
eS[eS>=1] = 0.99
S <- sapply(eS,rbinom,n = 1, size = 1)
S.ind <- which(S==1)

##--------------------------------------------------------------------
## RCT data
n.t <- length(S.ind)
X.t <- X[S.ind,]
A.t <- rbinom(n.t,1,0.5) # randomly assign treatment to trial participant

## Real world data
X <- matrix(rnorm(n2*p,1,0.5),ncol = p)
X <- scale(X) + 1
P.ind <- sample(1:n2,size = n.p) ## RWD id
X.p <- X[P.ind,]
eA <- trt_prop(X.p) # treatment propensity score
A.p <- sapply(eA,rbinom,n=1,size=1)

# Outcomes
Y.t <- outcome_mod_lognorm(X.t,A.t)
Y.p <- outcome_mod_lognorm(X.p,A.p)


eta.vec <- c(0,0.005,0.01,0.02,0.03,0.05)
B <- 50
fit.gaussian <- genRCT(X.t, X.p, A.t, A.p, Y.t, Y.p, family = 'gaussian', SBW.est = 'SBW-1', eta.vec = eta.vec, select.iter = 5, select.fix = TRUE, B = B, 
                       conf.level = 0.05, seed = NULL, boot.progress = 'number')
fit.gaussian


#####################
#  binary outcome  #
#####################

p <- 4
n.p <- 2000
n <- 2e4 ## population size
n2 <-2e5 ## population size

beta0 <- 5.3

## Generating data
set.seed(12345)
X <- matrix(rnorm(n*p,1,1),ncol = p)
X <- scale(X) + 1
eS <- sample_prop(X, beta0)

#print(sum(eS>=1))
eS[eS>=1] = 0.99
S <- sapply(eS,rbinom,n = 1, size = 1)
S.ind <- which(S==1)

##--------------------------------------------------------------------
## RCT data
n.t <- length(S.ind)
X.t <- X[S.ind,]
A.t <- rbinom(n.t,1,0.5) # randomly assign treatment to trial participant

## Real world data
X <- matrix(rnorm(n2*p,1,0.5),ncol = p)
X <- scale(X) + 1
P.ind <- sample(1:n2,size = n.p) ## RWD id
X.p <- X[P.ind,]
eA <- trt_prop(X.p) # treatment propensity score
A.p <- sapply(eA,rbinom,n=1,size=1)


# Outcomes
eY <- outcome_mod(X.t,A.t)
Y.t <- sapply(eY,rbinom,n = 1, size = 1)
eY <- outcome_mod(X.p,A.p)
Y.p <- sapply(eY,rbinom,n = 1, size = 1)

eta.vec <- c(0,0.005,0.01,0.02,0.03,0.05)
B <- 50
fit.binom <- genRCT(X.t, X.p, A.t, A.p, Y.t, Y.p, family = 'binomial', SBW.est = 'SBW-1', eta.vec = eta.vec, select.iter = 5, select.fix = TRUE, B = B, 
                    conf.level = 0.05, seed = NULL, boot.progress = 'number')
fit.binom
```

Estimators

  - Naive : the difference in sample means of the two treatment groups
    in the RCT sample

  - IPSW: the IPSW estimator

  - AIPSW: the AIPSW estimator

  - AIPSW(S) : the AIPSW estimator using the method of Sieves with
    \(g\)(X) = \(g_2\)(X) (SCAD method separately for \(\mu_a(X)\) and
    \(\pi_{\delta}(X)\)). The nuisance functions \(\mu_a(X), a=0, 1\)
    are estimated based on the trial sample

  - SBW-1 : SBW balancing means of all the original covariates

  - SBW-2 : SBW balancing means and decile indicators of the original
    covariates

  - SBW-3 : SBW balancing means of all original covariates and their
    squares and two-way interactions

SBW weights (Chattopadhyay et al. (2019)) for ATT (treatment : RWD,
control group : RCT) were first estimated using SBW-1, SBW-2, SBW-3
methods, and then plugged in to calculate IPSW.

  - CW: the CW estimator defined by with \(g\)(X) = \(g_1\)(X)

  - ACW-t: the ACW estimator defined by with \(g\)(X) = \(g_1\)(X) and
    the nuisance functions \(\mu_a(X), a=0, 1\) are estimated based on
    the trial sample

  - ACW-t(S): the penalized ACW-t estimator using the method of Sieves
    with \(g\)(X) = \(g_2\)(X) (“S” stands for method of sieves method)

  - ACW-t(S)-eff: the penalized ACW-t estimator to select sieve basis
    that are predictive of outcome and select those terms

  - ACW-b: the ACW estimator defined by with \(g\)(X) = \(g_1\)(X)and
    the nuisance functions \(\mu_a(X), a=0, 1\) are estimated are
    estimated based on both RCT and observational samples

  - ACW-b(S): the penalized ACW-b estimator using the method of sieves
    with \(g\)(X) = \(g_2\)(X)

  - ACW-t(S)-eff: the penalized ACW-b estimator to select sieve basis
    that are predictive of outcome and select those terms

where \(g_1(X) = (X_1, X_2, X_3, X_4)^T\) and
\(g_2(X) = (X_1, ... X_4, X_1X_2,...,X_3X_4,X_1^2,...,X_4^2)\)

``` 

knitr::kable(round(fit.gaussian$fit, 3), caption = 'ATE comparison for continuous outcomes') %>% kable_styling(bootstrap_options = "striped", full_width = T)
knitr::kable(round(fit.binom$fit, 3), caption = 'ATE comparison for binary outcomes') %>% kable_styling(bootstrap_options = "striped", full_width = T)
```
