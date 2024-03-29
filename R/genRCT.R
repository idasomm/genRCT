#' Generalizing the average treatment effect (ATE) from the trial using observational studies (OS)
#'
#' Provides estimation and inference for the average treatment effect (ATE), generalizing from a randomized controlled trial (RCT)
#' to a target population leveraging observational studies. The method of sieves can be used for the estimation of sampling score and outcome models.
#' The ATE can be estimated with up to 6 methods among the following: Naive, IPSW, AIPSW, CW, ACW-t, ACW-b.
#'
#' @param Y.trial Observed outcome from a trial; vector of size \code{n} (trial sample size).
#' @param X.trial Matrix of \code{p} baseline covariates from a trial; dimension \code{n} by \code{p}.
#' @param A.trial Treatment received from a trial; vector of size \code{n}.
#' @param Y.rwe Observed outcome from OS; if obtained, vector of size \code{m} (OS sample size);
#' otherwise, set \code{Y.rwe = NULL}.
#' @param X.rwe Matrix of \code{p} baseline covariates from OS; dimension \code{m} by \code{p}.
#' @param A.rwe Treatment received from OS; if obtained, vector of size \code{m}; otherwise, set \code{A.rwe = NULL}.
#' @param family The type of outcome; \code{'gaussian'} for gaussian regression or \code{'binomial'} for logistic regression.
#' Default is \code{"gaussian"}.
#' @param estimators A vector of one or multiple methods to estimate the ATE. Allowed values are
#' \code{'Naive'}, \code{'IPSW'}, \code{'AIPSW'}, \code{'CW'}, \code{'ACW-t'}, \code{'ACW-b'}.
#' The \code{'ACW-b'} is allowed only when both \code{Y.rwe} and \code{A.rwe} are obtained.
#' Default specifies all 6 methods.
#' @param sieve A logical value indicating whether the method of sieves are used for estimating sampling score and outcome models.
#' Used only if \code{estimators = 'AIPSW'} or \code{'ACW-t'} or \code{'ACW-b'}. Default is \code{TRUE}.
#' @param inference A logical value indicating whether inference for the ATE via bootstrap should be provided.
#' Default it \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{100}.
#' @param conf.level The level of bootstrap confidence interval; Default is \code{0.05}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param plot.boot A logical value indicating whether histograms of the bootstrap samples should be produced. Default is \code{TRUE}.
#' @param verbose A logical value indicating whether intermediate progress messages should be printed.
#' Default is \code{TRUE}.
#' @return Return a list containing:\tabular{ll}{
#'    \code{fit} \tab\tab  A table of estimated ATEs with bootstrap SE and confidence interval. \cr
#'    \tab \cr
#'    \code{plot} \tab\tab A set of histograms displaying the distribution of the bootstrapped estimates. The red vertical reference lines
#'     represent the estimated ATEs from each method. \cr
#' }
#'
#' @examples
#' fit <- genRCT(Y.trial = Y.trial, A.trial = A.trial, X.trial = X.trial, Y.rwe = Y.rwe, A.rwe = A.rwe, X.rwe = X.rwe,
#' family = "gaussian", estimators = c("Naive", "IPSW", "AIPSW", "CW", "ACW-t", "ACW-b"), sieve = TRUE,
#' inference = TRUE, n.boot = 500, conf.level = 0.05, seed = 123, plot.boot = TRUE, verbose = TRUE)
#'
#' @export
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom stats quantile model.matrix poly lm glm
#' @importFrom data.table data.table
#' @importFrom ncvreg cv.ncvreg
#' @importFrom nleqslv searchZeros

'genRCT' <- function(Y.trial, A.trial, X.trial, Y.rwe, A.rwe, X.rwe, family = "gaussian",
                     estimators = c("Naive", "IPSW", "AIPSW", "CW", "ACW-t", "ACW-b"), sieve = TRUE,
                     inference = TRUE, n.boot = 100, conf.level = 0.05, seed = NULL, plot.boot = TRUE, verbose = TRUE) {

  t.start <- Sys.time()

  # Check arguments
  if (is.null(A.rwe) | is.null(Y.rwe)) {
    if (c("ACW-b") %in% estimators) stop("The ACW-b estimator must have A.rwe and Y.rwe.")
  }

  # Change data.frame to matrix
  Y.trial <- c(as.matrix(Y.trial))
  A.trial <- c(as.matrix(A.trial))
  X.trial <- as.matrix(X.trial)
  X.rwe <- as.matrix(X.rwe)
  if (!is.null(Y.rwe)) Y.rwe <- c(as.matrix(Y.rwe))
  if (!is.null(A.rwe)) A.rwe <- c(as.matrix(A.rwe))

  # Fit estimators
  cat(" Fitting estimators.. \n")
  fit <- genRCT.estimators(Y.trial = Y.trial, A.trial = A.trial, X.trial = X.trial, Y.rwe = Y.rwe, A.rwe = A.rwe, X.rwe = X.rwe,
                           estimators = estimators, sieve = sieve, family = family, seed = seed,
                           osel1.t = NULL, osel0.t = NULL, osel1.b = NULL, osel0.b = NULL, osel.ipsw = NULL, flag.boot = FALSE)
  names(fit$ate) <- estimators

  if (inference == FALSE) {
    return(fit$ate)
  } else {
    n.est <- length(estimators)
    n1 <- nrow(X.trial)
    n2 <- nrow(X.rwe)
    tau_B <- matrix(NA, nrow = n.boot, ncol = n.est)
    if(is.null(seed)) seed <- 0

    cat(" Bootstrapping.. \n")
    if(verbose == TRUE) {
      pb <- txtProgressBar(min = 0, max = n.boot, style = 3)
    }
    for (i in 1:n.boot){
      samp1b <- sample(n1, replace = TRUE)
      samp2b <- sample(n2, replace = TRUE)
      X1b <- X.trial[samp1b, ]
      X2b <- X.rwe[samp2b, ]
      A1b <- A.trial[samp1b]
      A2b <- A.rwe[samp2b]
      Y1b <- Y.trial[samp1b]
      Y2b <- Y.rwe[samp2b]
      fit.boot <- genRCT.estimators(Y.trial = Y1b, A.trial = A1b, X.trial = X1b, Y.rwe = Y2b, A.rwe = A2b, X.rwe = X2b,
                                    estimators = estimators, sieve = sieve, family = family, seed = seed + i,
                                    osel1.t = fit$hypers$osel1.t, osel0.t = fit$hypers$osel0.t, osel1.b = fit$hypers$osel1.b,
                                    osel0.b = fit$hypers$osel0.b, osel.ipsw = fit$hypers$osel.ipsw, flag.boot = TRUE)
      tau_B[i,] <- unlist(fit.boot$ate)

      if (verbose == TRUE) {
        Sys.sleep(0.025)
        setTxtProgressBar(pb, i)
      }
    } # end of for(i in 1:n.boot)
    if (verbose == TRUE) close(pb)
    se <- apply(tau_B, 2, sd, na.rm = TRUE)
    ci <- apply(tau_B, 2, function(x) quantile(x, probs = c(conf.level / 2, 1 - conf.level / 2), na.rm = TRUE))

    est <- cbind(t(fit$ate), se, t(ci))
    colnames(est) <- c("ATE", "SE", colnames(est)[3:4])

    result <- c()
    result$fit <- est

    if(plot.boot == TRUE) {
      ate.df <- data.frame(fit$ate)
      colnames(ate.df) <- estimators
      ate.df <- ate.df %>% gather(key = "Estimators", value = "ATE") %>% mutate(Estimators = factor(Estimators, levels = estimators))
      df <- data.frame(tau_B)
      colnames(df) <- estimators
      p1 <- df %>% gather(key = "Estimators", value = "ATE") %>% na.omit() %>% mutate(Estimators = factor(Estimators, levels = estimators)) %>%
              ggplot(aes(x = ATE)) + geom_histogram(bins = 30, alpha = 0.7) + facet_wrap(. ~ Estimators, scales = "free") +
              geom_vline(data = ate.df, aes(xintercept = ATE), color = "red") +
              labs(x = "Bootstrap values", y = "Frequency", title = paste0(n.boot, " bootstrap iterations")) +
              theme_bw()
      print(p1)
      result$plot <- p1
    }

    t.end <- Sys.time()
    t.diff <- difftime(t.end, t.start)
    cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')

    return(result)
  } # end of if(inference == FALSE) {} else {}
}


