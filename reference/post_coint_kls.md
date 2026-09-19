# Posterior Draw for Cointegration Models

Produces a draw of coefficients for cointegration models with a prior on
the cointegration space as proposed in Koop et al. (2010) and a draw of
non-cointegration coefficients from a normal density.

## Usage

``` r
post_coint_kls(
  y,
  beta,
  w,
  sigma_i,
  v_i,
  p_tau_i,
  g_i,
  x = NULL,
  gamma_mu_prior = NULL,
  gamma_v_i_prior = NULL
)
```

## Arguments

- y:

  a \\K \times T\\ matrix of differenced endogenous variables.

- beta:

  a \\M \times r\\ cointegration matrix \\\beta\\.

- w:

  a \\M \times T\\ matrix of variables in the cointegration term.

- sigma_i:

  an inverse of the \\K \times K\\ variance-covariance matrix.

- v_i:

  a numeric between 0 and 1 specifying the shrinkage of the
  cointegration space prior.

- p_tau_i:

  an inverted \\M \times M\\ matrix specifying the central location of
  the cointegration space prior of \\sp(\beta)\\.

- g_i:

  a \\K \times K\\ matrix.

- x:

  a \\N \times T\\ matrix of differenced regressors and unrestricted
  deterministic terms.

- gamma_mu_prior:

  a \\KN \times 1\\ prior mean vector of non-cointegration coefficients.

- gamma_v_i_prior:

  an inverted \\KN \times KN\\ prior covariance matrix of
  non-cointegration coefficients.

## Value

A named list containing the following elements:

- alpha:

  a draw of the \\K \times r\\ loading matrix.

- beta:

  a draw of the \\M \times r\\ cointegration matrix.

- Pi:

  a draw of the \\K \times M\\ cointegration matrix \\\Pi = \alpha
  \beta^{\prime}\\.

- Gamma:

  a draw of the \\K \times N\\ coefficient matrix for non-cointegration
  parameters.

## Details

The function produces posterior draws of the coefficient matrices
\\\alpha\\, \\\beta\\ and \\\Gamma\\ for the model \$\$y\_{t} = \alpha
\beta^{\prime} w\_{t-1} + \Gamma z\_{t} + u\_{t},\$\$ where \\y\_{t}\\
is a K-dimensional vector of differenced endogenous variables.
\\w\_{t}\\ is an \\M \times 1\\ vector of variables in the cointegration
term, which include lagged values of endogenous and exogenous variables
in levels and restricted deterministic terms. \\z\_{t}\\ is an
N-dimensional vector of differenced endogenous and exogenous explanatory
variabes as well as unrestricted deterministic terms. The error term is
\\u_t \sim \Sigma\\.

Draws of the loading matrix \\\alpha\\ are obtained using the prior on
the cointegration space as proposed in Koop et al. (2010). The posterior
covariance matrix is \$\$\overline{V}\_{\alpha} = \left\[\left(v^{-1}
(\beta^{\prime} P\_{\tau}^{-1} \beta) \otimes G\_{-1}\right) +
\left(ZZ^{\prime} \otimes \Sigma^{-1} \right) \right\]^{-1}\$\$ and the
posterior mean by \$\$\overline{\alpha} = \overline{V}\_{\alpha} +
vec(\Sigma^{-1} Y Z^{\prime}),\$\$ where \\Y\\ is a \\K \times T\\
matrix of differenced endogenous variables and \\Z = \beta^{\prime} W\\
with \\W\\ as an \\M \times T\\ matrix of variables in the cointegration
term.

For a given prior mean vector \\\underline{\Gamma}\\ and prior
covariance matrix \\\underline{V\_{\Gamma}}\\ the posterior covariance
matrix of non-cointegration coefficients in \\\Gamma\\ is obtained by
\$\$\overline{V}\_{\Gamma} = \left\[ \underline{V}\_{\Gamma}^{-1} +
\left(X X^{\prime} \otimes \Sigma^{-1} \right) \right\]^{-1}\$\$ and the
posterior mean by \$\$\overline{\Gamma} = \overline{V}\_{\Gamma} \left\[
\underline{V}\_{\Gamma}^{-1} \underline{\Gamma} + vec(\Sigma^{-1} Y
X^{\prime}) \right\],\$\$ where \\X\\ is an \\M \times T\\ matrix of
explanatory variables, which do not enter the cointegration term.

Draws of the cointegration matrix \\\beta\\ are obtained using the prior
on the cointegration space as proposed in Koop et al. (2010). The
posterior covariance matrix of the unrestricted cointegration matrix
\\B\\ is \$\$\overline{V}\_{B} = \left\[\left(A^{\prime} G^{-1} A
\otimes v^{-1} P\_{\tau}^{-1} \right) + \left(A^{\prime} \Sigma^{-1} A
\otimes WW^{\prime} \right) \right\]^{-1}\$\$ and the posterior mean by
\$\$\overline{B} = \overline{V}\_{B} + vec(W Y\_{B}^{-1} \Sigma^{-1}
A),\$\$ where \\Y\_{B} = Y - \Gamma X\\ and \\A = \alpha
(\alpha^{\prime} \alpha)^{-\frac{1}{2}}\\.

The final draws of \\\alpha\\ and \\\beta\\ are calculated using \\\beta
= B (B^{\prime} B)^{-\frac{1}{2}}\\ and \\\alpha = A (B^{\prime}
B)^{\frac{1}{2}}\\.

## References

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224-242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

## Examples

``` r
 
# Load data
data("e6")

# Generate model data
temp <- create_bvecmodel(e6, p = 2, r = 1)
y <- t(temp$data$train$y)
ect <- t(temp$data$train$w)
x <- t(temp$data$train$x)

k <- nrow(y) # Endogenous variables
tt <- ncol(y) # Number of observations
r <- temp$model$rank # Cointegration rank

# Priors. They cover the alpha block as well as the Gamma coefficients,
# so the length is k * r + k * nrow(x).
n_ag <- k * (r + nrow(x))
gamma_mu_prior <- matrix(0, n_ag)
gamma_v_i_prior <- diag(0, n_ag)

# Initial value of Sigma
sigma <- tcrossprod(y) / tt
sigma_i <- solve(sigma)

# Initial values of beta
beta <- matrix(c(1, -4), nrow(ect))

# Draw parameters
coint <- post_coint_kls(y = y, beta = beta, w = ect, x = x, sigma_i = sigma_i,
                        v_i = 0, p_tau_i = diag(1, nrow(ect)), g_i = sigma_i,
                        gamma_mu_prior = gamma_mu_prior,
                        gamma_v_i_prior = gamma_v_i_prior)
```
