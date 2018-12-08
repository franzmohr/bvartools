
bvartools
=========

[![Build Status](https://travis-ci.org/franzmohr/bvartools.svg?branch=master)](https://travis-ci.org/franzmohr/bvartools)

Overview
--------

The package `bvartools` implements some common functions used for Bayesian inference for mulitvariate time series models. The package should give researchers maximum freedom in setting up a Gibbs sampler in R and keep calculation time limited at the same time. This is achieved by implementing posterior simultion functions in C++. Its main features are

-   Posterior simulation functions, which are written in C++ for faster calculation
-   The `bvars` function collects the output of a Gibbs sampler in a standardised object, which can be used for further analyses
-   Functions for further analyses such as `irf` for impulse response analysis

Installation
------------

### Development version

``` r
# install.packages("devtools")
devtools::install_github("franzmohr/bvartools")
```

Usage
-----

This example covers the estimation of a simple BVAR model. For further examples on time varying parameter (TVP), stochastic volatility (SV), and vector error correction models as well as shrinkage methods like stochastic search variable selection (SSVS) or Bayesian variable selection (BVS) see the vignettes of the package.

### Generate an artificial sample

The artifical data set for this example is based on the process

$$y\_{t} = \\begin{pmatrix} 0.4 & 0.3 & 0 \\\\ 0 & 0.6 & 0.1 \\\\ 0 & 0.2 & 0.2 \\end{pmatrix} y\_{t-1} + u\_t, \\ \\ \\text{with} \\ \\ u\_t \\sim N(0, \\Sigma) \\ \\ \\text{and} \\ \\ \\Sigma = \\begin{pmatrix} 0.09 & 0 & 0 \\\\ 0 & 0.09 & 0 \\\\ 0 & 0 & 0.09 \\end{pmatrix}.$$

``` r
# Number of endogenous variables
k <- 3

# Define A
A <- matrix(c(.4, 0, 0, .3, .6, .2, 0, .1, .2), k)

# Number of observations
t <- 150

# Set seed of random number generator for reproducibility
set.seed(1234567)

# Generate matrix with empty values
y <- matrix(NA, 3, t + 1)

# Initial value of the series
y[, 1] <- rnorm(k, 0, .3)

# Recursively calculate values of the series
for (i in 2:(t + 1)) {
  y[, i] <- A %*% y[, i - 1] + rnorm(k, 0, .3)
}

data <- ts(t(y)) # Transform into a 'time series' object
dimnames(data)[[2]] <- c("s1", "s2", "s3") # Rename variables

plot(data) # Plot the series
```

<img src="README_files/figure-markdown_github/art-data-1.png" style="display: block; margin: auto;" />

### Prepare data for estimation

``` r
library(bvartools)

temp <- gen_var(data, p = 1, deterministic = "none")

y <- temp$y
x <- temp$x
```

### Estimation

``` r
iter <- 10000 # Number of iterations of the Gibbs sampler
burnin <- 2000 # Number of burn-in draws

t <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
nvars <- k * nrow(x) # Number of estimated coefficients

# Set priors
A_mu_prior <- matrix(0, nvars) # Vector of prior parameter means
A_V_i_prior <- diag(0, nvars) # Inverse of the prior covariance matrix

Sigma_df_prior <- k # Prior degrees of freedom
Sigma_V_prior <- diag(.00001, k) # Prior covariance matrix
Sigma_df_post <- t + Sigma_df_prior # Posterior degrees of freedom

# Initial values
Sigma_i_draw <- rWishart(1, k, solve(Sigma_V_prior))[,,1]
Sigma_draw <- solve(Sigma_i_draw)

# Data containers
store <- iter - burnin
draws_A <- matrix(NA, nvars, store)
draws_Sigma <- matrix(NA, k^2, store)
draws_LL <- matrix(NA, t, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  A_draw <- post_normal(y, x, Sigma_i_draw, A_mu_prior, A_V_i_prior)
  
  # Draw variance-covariance matrix
  res <- y - matrix(A_draw, k) %*% x # Obtain residuals
  Sigma_V_post <- solve(Sigma_V_prior + tcrossprod(res))
  Sigma_i_draw <- rWishart(1, Sigma_df_post, Sigma_V_post)[,, 1]
  Sigma_draw <- solve(Sigma_i_draw) # Invert Sigma_i to obtain Sigma
  
  # Store draws
  if (draw > burnin) {
    draws_A[, draw - burnin] <- A_draw
    draws_Sigma[, draw - burnin] <- Sigma_draw
    
    # Calculate Log-Likelihood
    draws_LL[, draw - burnin] <- loglik_gauss(res, Sigma_draw, Sigma_i_draw)
  }
}
```

Obtain point estimates as the mean of the parameter draws

``` r
A <- rowMeans(draws_A) # Obtain means for every row
A <- matrix(A, k) # Transform mean vector into a matrix
A <- round(A, 2) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

A # Print
```

    ##     s1.1 s2.1  s3.1
    ## s1  0.46 0.28 -0.09
    ## s2 -0.04 0.70  0.14
    ## s3  0.06 0.24  0.17

``` r
Sigma <- rowMeans(draws_Sigma) # Obtain means for every row
Sigma <- matrix(Sigma, k) # Transform mean vector into a matrix
Sigma <- round(Sigma, 2) # Round values
dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions

Sigma # Print
```

    ##       s1  s2    s3
    ## s1  0.08 0.0 -0.01
    ## s2  0.00 0.1  0.00
    ## s3 -0.01 0.0  0.09

The means of the coefficient draws are very close to the results of the frequentist estimatior and, hence, also close to the true parameter values.

### `bvars` objects

The `bvars` function can be used to collect relevant output of the Gibbs sampler into a standardised object, which can be used by further function such as `irf` to obtain impulse responses.

``` r
bvar_est <- bvars(y = y, x = x, A = draws_A, Sigma = draws_Sigma)
```

### Impulse response analysis

#### Forecast error impulse response

``` r
IR <- irf(bvar_est, impulse = "s2", response = "s1", n.ahead = 15)

plot(IR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-markdown_github/feir-1.png)

#### Orthogonalised impulse response

``` r
OIR <- irf(bvar_est, impulse = "s2", response = "s1", n.ahead = 15, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-markdown_github/oir-1.png)

#### Generalised impulse response

``` r
GIR <- irf(bvar_est, impulse = "s2", response = "s1", n.ahead = 15, type = "gir")

plot(GIR, main = "Generalised Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-markdown_github/gir-1.png) \#\# References

Lütkepohl, H. (2007). *New introduction to multiple time series analyis*. Berlin: Springer.

Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. *Economics Letters*, 58, 17-29.
