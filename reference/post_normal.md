# Posterior Draw from a Normal Distribution

Produces a draw of coefficients from a normal posterior density.

## Usage

``` r
post_normal(y, x, sigma_i, a_prior, v_i_prior)
```

## Arguments

- y:

  a \\K \times T\\ matrix of endogenous variables.

- x:

  an \\M \times T\\ matrix of explanatory variables.

- sigma_i:

  the inverse of the \\K \times K\\ variance-covariance matrix.

- a_prior:

  a \\KM \times 1\\ numeric vector of prior means.

- v_i_prior:

  the inverse of the \\KM \times KM\\ prior covariance matrix.

## Value

A vector.

## Details

The function produces a vectorised posterior draw \\a\\ of the \\K
\times M\\ coefficient matrix \\A\\ for the model \$\$y\_{t} = A
x\_{t} + u\_{t},\$\$ where \\y\_{t}\\ is a K-dimensional vector of
endogenous variables, \\x\_{t}\\ is an M-dimensional vector of
explanatory variabes and the error term is \\u_t \sim \Sigma\\.

For a given prior mean vector \\\underline{a}\\ and prior covariance
matrix \\\underline{V}\\ the posterior covariance matrix is obtained by
\$\$\overline{V} = \left\[ \underline{V}^{-1} + \left(X X^{\prime}
\otimes \Sigma^{-1} \right) \right\]^{-1}\$\$ and the posterior mean by
\$\$\overline{a} = \overline{V} \left\[ \underline{V}^{-1}
\underline{a} + vec(\Sigma^{-1} Y X^{\prime}) \right\],\$\$ where \\Y\\
is a \\K \times T\\ matrix of the endogenous variables and \\X\\ is an
\\M \times T\\ matrix of the explanatory variables.

## References

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## Examples

``` r

# Load data
data("e1")
data <- diff(log(e1))

# Generate model data
temp <- create_bvarmodel(data, p = 2, deterministic = "const")
y <- t(temp$data$train$y)
x <- t(temp$data$train$x)
k <- nrow(y)
tt <- ncol(y)
m <- k * nrow(x)

# Priors
a_mu_prior <- matrix(0, m)
a_v_i_prior <- diag(0.1, m)

# Initial value of inverse Sigma
sigma_i <- solve(tcrossprod(y) / tt)

# Draw parameters
a <- post_normal(y = y, x = x, sigma_i = sigma_i,
                 a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
```
