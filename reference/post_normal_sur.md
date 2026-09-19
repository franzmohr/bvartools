# Posterior Draw from a Normal Distribution

Produces a draw of coefficients from a normal posterior density for a
model with seemingly unrelated regresssions (SUR).

## Usage

``` r
post_normal_sur(y, z, sigma_i, a_prior, v_i_prior, svd = FALSE)
```

## Arguments

- y:

  a \\K \times T\\ matrix of endogenous variables.

- z:

  a \\KT \times M\\ matrix of explanatory variables.

- sigma_i:

  the inverse of the constant \\K \times K\\ error variance-covariance
  matrix. For time varying variance-covariance matrics a \\KT \times K\\
  can be provided.

- a_prior:

  a \\M x 1\\ numeric vector of prior means.

- v_i_prior:

  the inverse of the \\M x M\\ prior covariance matrix.

- svd:

  logical. If `TRUE` the singular value decomposition is used to
  determine the root of the posterior covariance matrix. Default is
  `FALSE` which means that the eigenvalue decomposition is used.

## Value

A vector.

## Details

The function produces a posterior draw of the coefficient vector \\a\\
for the model \$\$y\_{t} = Z\_{t} a + u\_{t},\$\$ where \\u_t \sim N(0,
\Sigma\_{t})\\. \\y_t\\ is a K-dimensional vector of endogenous
variables and \\Z_t = z_t^{\prime} \otimes I_K\\ is a \\K \times KM\\
matrix of regressors with \\z_t\\ as a vector of regressors.

For a given prior mean vector \\\underline{a}\\ and prior covariance
matrix \\\underline{V}\\ the posterior covariance matrix is obtained by
\$\$\overline{V} = \left\[ \underline{V}^{-1} + \sum\_{t=1}^{T}
Z\_{t}^{\prime} \Sigma\_{t}^{-1} Z\_{t} \right\]^{-1}\$\$ and the
posterior mean by \$\$\overline{a} = \overline{V} \left\[
\underline{V}^{-1} \underline{a} + \sum\_{t=1}^{T} Z\_{t}^{\prime}
\Sigma\_{t}^{-1} y\_{t} \right\].\$\$

## Examples

``` r

# Load data
data("e1")
data <- diff(log(e1))

# Generate model data
temp <- create_bvarmodel(data, p = 2, deterministic = "const")
y <- t(temp$data$train$y)
z <- temp$data$train$z
k <- nrow(y)
tt <- ncol(y)
m <- ncol(z)

# Priors
a_mu_prior <- matrix(0, m)
a_v_i_prior <- diag(0.1, m)

# Initial value of inverse Sigma
sigma_i <- solve(tcrossprod(y) / tt)

# Draw parameters
a <- post_normal_sur(y = y, z = z, sigma_i = sigma_i,
                     a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
```
