# Stochastic Search Variable Selection

`ssvs` employs stochastic search variable selection as proposed by
George et al. (2008) to produce a draw of the precision matrix of the
coefficients in a VAR model.

## Usage

``` r
ssvs(a, tau0, tau1, prob_prior, include = NULL)
```

## Arguments

- a:

  an M-dimensional vector of coefficient draws.

- tau0:

  an M-dimensional vector of prior standard deviations for restricted
  coefficients in vector `a`.

- tau1:

  an M-dimensional vector of prior standard deviations for unrestricted
  coefficients in vector `a`.

- prob_prior:

  an M-dimensional vector of prior inclusion probabilites for the
  coefficients in vector `a`.

- include:

  an integer vector specifying the positions of coefficients in vector
  `a`, which should be included in the SSVS algorithm. If `NULL`
  (default), SSVS will be applied to all coefficients.

## Value

A named list containing two components:

- v_i:

  an \\M \times M\\ inverse prior covariance matrix.

- lambda:

  an M-dimensional vector of inclusion parameters.

## Details

The function employs stochastic search variable selection (SSVS) as
proposed by George et al. (2008) to produce a draw of the diagonal
inverse prior covariance matrix \\\underline{V}^{-1}\\ and the
corresponding vector of inclusion parameters \\\lambda\\ of the
vectorised coefficient matrix \\a = vec(A)\\ for the VAR model \$\$y_t =
A x_t + u_t,\$\$ where \\y\_{t}\\ is a K-dimensional vector of
endogenous variables, \\x\_{t}\\ is a vector of explanatory variabes and
the error term is \\u_t \sim \Sigma\\.

## References

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

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

# Obtain SSVS priors using the semiautomatic approach
priors <- ssvs_prior(temp, semiautomatic = c(0.1, 10))
tau0 <- priors$tau0
tau1 <- priors$tau1

# Prior for inclusion parameter
prob_prior <- matrix(0.5, m)

# Priors
a_mu_prior <- matrix(0, m)
a_v_i_prior <- diag(c(tau1^2), m)

# Initial value of Sigma
sigma_i <- solve(tcrossprod(y) / tt)

# Draw parameters
a <- post_normal(y = y, x = x, sigma_i = sigma_i,
                 a_prior = a_mu_prior, v_i_prior = a_v_i_prior)

# Run SSVS
lambda <- ssvs(a = a, tau0 = tau0, tau1 = tau1,
               prob_prior = prob_prior)
```
