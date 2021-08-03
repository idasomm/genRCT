
# genRCT

### Generalizing the average treatment effect (ATE) from the trial using observational studies (OS)

### Installation

``` r
devtools::install_github("idasomm/genRCT", 
                         auth_token = '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')
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
| Y.trial    | Observed outcome from a trial; vector of size n (trial sample size).                                                                                                                                                                   |
| A.trial    | Treatment received from a trial; vector of size n.                                                                                                                                                                                     |
| X.trial    | Matrix of p baseline covariates from a trial; dimension n by p.                                                                                                                                                                        |
| Y.rwe      | Observed outcome from OS; if obtained, vector of size m (OS sample size); otherwise, set Y.rwe = NULL.                                                                                                                                 |
| A.rwe      | Treatment received from OS; if obtained, vector of size m; otherwise, set A.rwe = NULL.                                                                                                                                                |
| X.rwe      | Matrix of p baseline covariates from OS; dimension m by p.                                                                                                                                                                             |
| family     | The type of outcome; “gaussian” for gaussian regression or “binomial” for logistic regression Default is “gaussian”.                                                                                                                   |
| estimators | A vector of one or multiple methods to estimate the ATE. Allowed values are “Naive”, “IPSW”, “AIPSW”, “CW”, “ACW-t”, “ACW-b”. The “ACW-b” is allowed only when both “Y.rwe” and “A.rwe” are obtained. Default specifies all 6 methods. |
| sieve      | A logical value indicating whether the method of sieves are used for estimating sampling score and outcome models. Used only if estimators = “AIPSW” or “ACW-t” or “ACW-b”. Default is TRUE.                                           |
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

``` r

library(genRCT)
library(tidyverse)
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
#> ✓ tibble  3.1.0     ✓ dplyr   1.0.5
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following objects are masked from 'package:dplyr':
#> 
#>     between, first, last
#> The following object is masked from 'package:purrr':
#> 
#>     transpose

data("simulData")

data.trial <- simulData %>% dplyr::filter(delta == 1)
data.rwe <- simulData %>% dplyr::filter(delta == 0)

Y.trial <- data.trial %>% dplyr::select(Y)
A.trial <- data.trial %>% dplyr::select(A)
X.trial <- data.trial %>% dplyr::select(X1:X5)

Y.rwe <- data.rwe %>% dplyr::select(Y)
A.rwe <- data.rwe %>% dplyr::select(A)
X.rwe <- data.rwe %>% dplyr::select(X1:X5)

fit <- genRCT(Y.trial = Y.trial, A.trial = A.trial, X.trial = X.trial, Y.rwe = Y.rwe, 
              A.rwe = A.rwe, X.rwe = X.rwe, family = "gaussian",
              estimators = c("Naive", "IPSW", "AIPSW", "CW", "ACW-t", "ACW-b"), 
              sieve = TRUE, inference = TRUE, n.boot = 500, conf.level = 0.05, 
              seed = 12345, plot.boot = TRUE, verbose = FALSE)
#>  Fitting estimators.. 
#>  Bootstrapping..
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

    #> Total runtime : 5.622568 mins
    
    fit$fit
    #>            ATE         SE      2.5%    97.5%
    #> Naive 50.31553  2.4868271 45.342828 55.36317
    #> IPSW  35.69669  8.1734416 19.479857 49.27143
    #> AIPSW 26.63354  0.6464247 25.283526 27.78797
    #> CW    27.74951 14.7089762  1.937296 52.01549
    #> ACW-t 26.88160  0.7517235 25.197760 28.06301
    #> ACW-b 27.57737  0.7093065 25.914611 28.65249
    
    truth <- 27.4
    fit$plot + geom_vline(xintercept = truth, color = 'blue', linetype = "dashed")

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" /> The
red solid lines represent the estimated ATEs and the blue dashed lines
represent the true ATE. The differences are

    #>       Naive        IPSW       AIPSW          CW       ACW-t       ACW-b 
    #> -22.9155252  -8.2966905   0.7664578  -0.3495104   0.5184023  -0.1773662
