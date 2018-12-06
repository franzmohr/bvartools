
bvartools
=========

[![Build Status](https://travis-ci.org/franzmohr/bvartools.svg?branch=master)](https://travis-ci.org/franzmohr/bvartools)

Overview
--------

bvartools implements some standard functions for Bayesian inference for mulitvariate time series analysis.

Installation
------------

### Development version

``` r
# install.packages("devtools")
devtools::install_github("franzmohr/bvartools")
```

Usage
-----

### Generate an artificial sample

``` r
# Define A
A <- matrix(c(.4, 0, 0, .3, .6, .2, 0, .1, .2), 3)

# Number of observations
t <- 150

# Set seed of random number generator for reproducibility
set.seed(1234567)

# Generate matrix with empty values
y <- matrix(NA, 3, t + 1)

# Initial value of the series
y[, 1] <- rnorm(3, 0, .3)

# Recursively calculate values of the series
for (i in 2:(t + 1)) {
  y[, i] <- A %*% y[, i - 1] + rnorm(3, 0, .3)
}

data <- ts(t(y)) # Transform into a 'time series' object
dimnames(data)[[2]] <- c("s1", "s2", "s3") # Rename variables

plot(data) # Plot the series
```

<img src="README_files/figure-markdown_github/art-data-1.png" style="display: block; margin: auto;" />

### Prepare data for estimation

``` r
p <- 1 # Number of lags
raw <- data
raw_names <- data_names <- dimnames(raw)[[2]]
for (i in 1:p) {
  data <- cbind(data, lag(raw, -i))
  data_names <- c(data_names, paste(raw_names, ".l", i, sep = ""))
}
data <- na.omit(data)
dimnames(data)[[2]] <- data_names

y <- t(data[, 1:3])
x <- t(data[, -(1:3)])
```

We will only require variables `y` and `x` in the following.

### Estimation

``` r
# Load the package
library(bvartools)

iter <- 10000 # Number of iterations of the Gibbs sampler
burnin <- 2000 # Number of burn-in draws

t <- ncol(y) # Number of observations
n <- nrow(y) # Number of endogenous variables
nvars <- n * nrow(x) # Number of estimated coefficients

# Set priors
A_mu_prior <- matrix(0, nvars) # Vector of prior parameter means
A_V_i_prior <- diag(0, nvars) # Inverse of the prior covariance matrix

Sigma_df_prior <- n # Prior degrees of freedom
Sigma_V_prior <- diag(.00001, n) # Prior covariance matrix
Sigma_df_post <- t + Sigma_df_prior # Posterior degrees of freedom

# Initial values
Sigma_i_draw <- rWishart(1, n, solve(Sigma_V_prior))[,,1]
Sigma_draw <- solve(Sigma_i_draw)

# Data containers
store <- iter - burnin
draws_A <- matrix(NA, nvars, store)
draws_Sigma <- matrix(NA, n^2, store)
draws_LL <- matrix(NA, t, store)

# Start Gibbs sampler
for (draw in 1:iter) {
  # Draw conditional mean parameters
  A_draw <- post_normal(y, x, Sigma_i_draw, A_mu_prior, A_V_i_prior)
  
  # Draw variance-covariance matrix
  res <- y - matrix(A_draw, n) %*% x # Obtain residuals
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
A <- matrix(A, n) # Transform mean vector into a matrix
A <- round(A, 2) # Round values
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions

A # Print
```

    ##    s1.l1 s2.l1 s3.l1
    ## s1  0.46  0.28 -0.09
    ## s2 -0.04  0.70  0.14
    ## s3  0.06  0.24  0.17

``` r
Sigma <- rowMeans(draws_Sigma) # Obtain means for every row
Sigma <- matrix(Sigma, n) # Transform mean vector into a matrix
Sigma <- round(Sigma, 2) # Round values
dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions

Sigma # Print
```

    ##       s1  s2    s3
    ## s1  0.08 0.0 -0.01
    ## s2  0.00 0.1  0.00
    ## s3 -0.01 0.0  0.09

The means of the coefficient draws are very close to the results of the frequentist estimatior and, hence, also close to the true parameter values.

### Impulse response

``` r
temp <- bvars(y = y, x = x, A = draws_A, Sigma = draws_Sigma)

IR <- irf(temp, impulse = "s2", response = "s1", n.ahead = 15)

plot(IR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
```

![](README_files/figure-markdown_github/ir-1.png)
