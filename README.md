
# bvartools

[![CRAN
status](https://www.r-pkg.org/badges/version/bvartools)](https://cran.r-project.org/package=bvartools)
[![Travis build
status](https://travis-ci.org/franzmohr/bvartools.svg?branch=master)](https://travis-ci.org/franzmohr/bvartools)

## Overview

The package `bvartools` implements some common functions used for
Bayesian inference for mulitvariate time series models. It should give
researchers maximum freedom in setting up an MCMC algorithm in R and
keep calculation time limited at the same time. This is achieved by
implementing posterior simulation functions in C++. Its main features
are

  - The `bvar` and `bvec` function collects the output of a Gibbs
    sampler in standardised objects, which can be used for further
    analyses
  - Further functions such as `predict`, `irf`, `fevd` for forecasting,
    impulse response analysis and forecast error variance decomposition,
    respectively.
  - Computationally intensive functions - such as for posterior
    simulation - are written in C++ using the `RcppArmadillo` package of
    Eddelbuettel and Sanderson (2014).\[1\]

Similar packages worth checking out are

  - [BVAR](https://cran.r-project.org/package=BVAR)
  - [bvarsv](https://cran.r-project.org/package=bvarsv)
  - [bvar](https://github.com/nk027/bvar)
  - [bvarr](https://github.com/bdemeshev/bvarr)
  - [bvars](https://github.com/joergrieger/bvars)
  - [mfbvar](https://github.com/ankargren/mfbvar)
  - [BMR](https://github.com/kthohr/BMR)

## Installation

``` r
install.packages("bvartools")
```

### Development version

``` r
# install.packages("devtools")
devtools::install_github("franzmohr/bvartools")
```

## Usage

This example covers the estimation of a simple Bayesian VAR (BVAR)
model. For further examples on time varying parameter (TVP), stochastic
volatility (SV), and vector error correction (VEC) models as well as
shrinkage methods like stochastic search variable selection (SSVS) or
Bayesian variable selection (BVS) see the vignettes of the package and
[r-econometrics.com](https://www.r-econometrics.com/timeseriesintro/).

### Data

To illustrate the estimation process the dataset E1 from Lütkepohl
(2007) is used. It contains data on West German fixed investment,
disposable income and consumption expenditures in billions of DM from
1960Q1 to 1982Q4.

``` r
library(bvartools)

data("e1")
e1 <- diff(log(e1))

plot(e1) # Plot the series
```

<img src="README_files/figure-gfm/data-1.png" style="display: block; margin: auto;" />

### Prepare data for estimation

The `gen_var` function produces the inputs `y` and `x` for the BVAR
estimator, where `y` is the matrix of dependent variables and `x` is the
matrix of regressors.

``` r
data <- gen_var(e1, p = 2, deterministic = "const")

y <- data$Y[, 1:73]
x <- data$Z[, 1:73]
```

As in Lütkepohl (2007) only the first 73 observations are used.

### Estimation

The following code sets up a simple Gibbs sampler algorithm.

``` r
set.seed(1234567)

iter <- 15000 # Number of iterations of the Gibbs sampler
burnin <- 5000 # Number of burn-in draws
store <- iter - burnin

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Set (uninformative) priors
a_mu_prior <- matrix(0, m) # Vector of prior parameter means
a_v_i_prior <- diag(0, m) # Inverse of the prior covariance matrix

u_sigma_df_prior <- 0 # Prior degrees of freedom
u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
u_sigma_i <- diag(.00001, k)
u_sigma <- solve(u_sigma_i)

# Data containers for posterior draws
draws_a <- matrix(NA, m, store)
draws_sigma <- matrix(NA, k^2, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  
  # Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x # Obtain residuals
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    draws_sigma[, draw - burnin] <- u_sigma
  }
}
```

### `bvar` objects

The `bvar` function can be used to collect relevant output of the Gibbs
sampler in a standardised object, which can be used by further functions
such as `predict` to obtain forecasts or `irf` for impulse respons
analysis.

``` r
bvar_est <- bvar(y = y, x = x, A = draws_a[1:18,],
                 C = draws_a[19:21, ], Sigma = draws_sigma)
```

``` r
summary(bvar_est)
```

    ## 
    ## Model:
    ## 
    ## y ~ invest.1 + income.1 + cons.1 + invest.2 + income.2 + cons.2 + const
    ## 
    ## Variable: invest 
    ## 
    ##             Mean      SD  Naive SD Time-series SD     2.5%      50%    97.5%
    ## invest.1 -0.3210 0.12910 0.0012910      0.0012963 -0.57576 -0.32234 -0.06704
    ## income.1  0.1468 0.56920 0.0056920      0.0056920 -0.97536  0.14389  1.27628
    ## cons.1    0.9661 0.68461 0.0068461      0.0066912 -0.37530  0.95599  2.32208
    ## invest.2 -0.1601 0.12744 0.0012744      0.0013199 -0.40851 -0.16121  0.09073
    ## income.2  0.1036 0.55673 0.0055673      0.0055673 -0.98186  0.09850  1.22167
    ## cons.2    0.9359 0.69796 0.0069796      0.0069975 -0.43081  0.92389  2.32856
    ## const    -0.0166 0.01774 0.0001774      0.0001774 -0.05164 -0.01661  0.01778
    ## 
    ## Variable: income 
    ## 
    ##               Mean       SD  Naive SD Time-series SD      2.5%       50%
    ## invest.1  0.043536 0.033070 3.307e-04      3.356e-04 -0.021750  0.043455
    ## income.1 -0.152419 0.143312 1.433e-03      1.443e-03 -0.435520 -0.151075
    ## cons.1    0.286935 0.173352 1.734e-03      1.734e-03 -0.055373  0.286012
    ## invest.2  0.049782 0.032382 3.238e-04      3.238e-04 -0.014160  0.049723
    ## income.2  0.018945 0.139841 1.398e-03      1.398e-03 -0.255458  0.017515
    ## cons.2   -0.008832 0.172217 1.722e-03      1.722e-03 -0.346456 -0.009778
    ## const     0.015782 0.004494 4.494e-05      4.569e-05  0.007125  0.015761
    ##            97.5%
    ## invest.1 0.10836
    ## income.1 0.12729
    ## cons.1   0.63148
    ## invest.2 0.11349
    ## income.2 0.29004
    ## cons.2   0.32871
    ## const    0.02462
    ## 
    ## Variable: cons 
    ## 
    ##               Mean       SD  Naive SD Time-series SD      2.5%       50%
    ## invest.1 -0.002661 0.026741 2.674e-04      2.674e-04 -0.055640 -0.002423
    ## income.1  0.223297 0.117474 1.175e-03      1.175e-03 -0.001396  0.222052
    ## cons.1   -0.262860 0.140987 1.410e-03      1.410e-03 -0.543567 -0.262607
    ## invest.2  0.033769 0.026286 2.629e-04      2.629e-04 -0.017495  0.034067
    ## income.2  0.354436 0.112782 1.128e-03      1.111e-03  0.130712  0.355690
    ## cons.2   -0.019983 0.139373 1.394e-03      1.361e-03 -0.291804 -0.019175
    ## const     0.012909 0.003664 3.664e-05      3.664e-05  0.005685  0.012894
    ##            97.5%
    ## invest.1 0.04963
    ## income.1 0.45050
    ## cons.1   0.01277
    ## invest.2 0.08634
    ## income.2 0.57461
    ## cons.2   0.25449
    ## const    0.02012
    ## 
    ## Variance-covariance matrix:
    ## 
    ##                    Mean        SD  Naive SD Time-series SD       2.5%       50%
    ## invest_invest 2.267e-03 4.118e-04 4.118e-06      4.615e-06  1.605e-03 2.215e-03
    ## invest_income 7.688e-05 7.493e-05 7.493e-07      8.259e-07 -6.562e-05 7.406e-05
    ## invest_cons   1.317e-04 6.171e-05 6.171e-07      6.961e-07  2.122e-05 1.291e-04
    ## income_invest 7.688e-05 7.493e-05 7.493e-07      8.259e-07 -6.562e-05 7.406e-05
    ## income_income 1.461e-04 2.673e-05 2.673e-07      2.928e-07  1.033e-04 1.428e-04
    ## income_cons   6.547e-05 1.732e-05 1.732e-07      1.923e-07  3.635e-05 6.381e-05
    ## cons_invest   1.317e-04 6.171e-05 6.171e-07      6.961e-07  2.122e-05 1.291e-04
    ## cons_income   6.547e-05 1.732e-05 1.732e-07      1.923e-07  3.635e-05 6.381e-05
    ## cons_cons     9.499e-05 1.723e-05 1.723e-07      1.922e-07  6.693e-05 9.285e-05
    ##                   97.5%
    ## invest_invest 0.0031925
    ## invest_income 0.0002346
    ## invest_cons   0.0002631
    ## income_invest 0.0002346
    ## income_income 0.0002078
    ## income_cons   0.0001053
    ## cons_invest   0.0002631
    ## cons_income   0.0001053
    ## cons_cons     0.0001338

The means of the posterior draws are very close to the results of the
frequentist estimatior in Lütkepohl (2007).

### Forecasts

Forecasts can be obtained with the function `predict`. If the model
contains deterministic terms, new values have to be provided in the
argument `new_D`, which must be of the same length as the argument
`n.ahead`.

``` r
bvar_pred <- predict(bvar_est, n.ahead = 10, new_D = rep(1, 10))

plot(bvar_pred)
```

![](README_files/figure-gfm/forecasts-1.png)<!-- -->

### Impulse response analysis

#### Forecast error impulse response

``` r
IR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8)

plot(IR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-gfm/feir-1.png)<!-- -->

#### Orthogonalised impulse response

``` r
OIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-gfm/oir-1.png)<!-- -->

#### Generalised impulse response

``` r
GIR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8, type = "gir")

plot(GIR, main = "Generalised Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-gfm/gir-1.png)<!-- -->

### Forecast error variance decomposition

``` r
bvar_fevd <- fevd(bvar_est, response = "cons")

plot(bvar_fevd, main = "FEVD of consumption")
```

![](README_files/figure-gfm/fevd-1.png)<!-- -->

## References

Eddelbuettel, D., & Sanderson C. (2014). RcppArmadillo: Accelerating R
with high-performance C++ linear algebra. *Computational Statistics and
Data Analysis, 71*, 1054-1063.
<https://doi.org/10.1016/j.csda.2013.02.005>

Lütkepohl, H. (2007). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis
in linear multivariate models. *Economics Letters, 58*, 17-29.
<https://doi.org/10.1016/S0165-1765(97)00214-0>

Sanderson, C., & Curtin, R. (2016). Armadillo: a template-based C++
library for linear algebra. *Journal of Open Source Software, 1*(2), 26.
<https://doi.org/10.21105/joss.00026>

1.  `RcppArmadillo` is the `Rcpp` bridge to the open source ‘Armadillo’
    library of Sanderson and Curtin (2016).
