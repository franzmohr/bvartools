# Bayesian Vector Autoregression Objects

`bvar` collects the posterior draws of a vector autoregressive model in
an object of class 'bvarmodel'.

## Usage

``` r
bvar(
  data = NULL,
  exogen = NULL,
  y,
  x = NULL,
  z = NULL,
  A0 = NULL,
  A = NULL,
  B = NULL,
  C = NULL,
  Sigma = NULL,
  error = NULL,
  varsel = NULL,
  iterations = NULL,
  burnin = 0
)
```

## Arguments

- data:

  the original time-series object of endogenous variables. If `NULL`
  (default), the object provided in argument `y` is used.

- exogen:

  the original time-series object of unmodelled variables. If `NULL`
  (default), it is reconstructed from the current values of those
  variables in argument `x`.

- y:

  a time-series object of endogenous variables with \\T\\ observations,
  usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- x:

  a time-series object of \\(pK + (1+s)M + N)\\ regressor variables,
  usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- z:

  a \\TK \times K (pK + (1+s)M + N)\\ data matrix, with \\K(K - 1)/2\\
  further columns for the contemporaneous endogenous variables of a
  structural model, usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).
  If `NULL` (default), it is generated from argument `x`.

- A0:

  either a \\K^2 \times S\\ matrix of MCMC coefficient draws of
  structural parameters or a named list, where element `coeffs` contains
  a \\K^2 \times S\\ matrix of MCMC coefficient draws of structural
  parameters and element `lambda` contains the corresponding draws of
  inclusion parameters in case variable selection algorithms were
  employed. For time varying parameter models the coefficient matrix
  must be \\TK^2 \times S\\. Draws of the error covariance matrix of the
  state equation can be provided as a \\K^2 \times S\\ matrix in an
  additional list element.

- A:

  either a \\pK^2 \times S\\ matrix of MCMC coefficient draws of lagged
  endogenous variables or a named list, where element `coeffs` contains
  a \\pK^2 \times S\\ matrix of MCMC coefficient draws of lagged
  endogenous variables and element `lambda` contains the corresponding
  draws of inclusion parameters in case variable selection algorithms
  were employed. For time varying parameter models the coefficient
  matrix must be \\pTK^2 \times S\\. Draws of the error covariance
  matrix of the state equation can be provided as a \\pK^2 \times S\\
  matrix in an additional list element.

- B:

  either a \\((1 + s)MK) \times S\\ matrix of MCMC coefficient draws of
  unmodelled, non-deterministic variables or a named list, where element
  `coeffs` contains a \\((1 + s)MK) \times S\\ matrix of MCMC
  coefficient draws of unmodelled, non-deterministic variables and
  element `lambda` contains the corresponding draws of inclusion
  parameters in case variable selection algorithms were employed. For
  time varying parameter models the coefficient matrix must be \\(1 +
  s)TMK \times S\\. Draws of the error covariance matrix of the state
  equation can be provided as a \\(1 + s)MK \times S\\ matrix in an
  additional list element.

- C:

  either a \\KN \times S\\ matrix of MCMC coefficient draws of
  deterministic terms or a named list, where element `coeffs` contains a
  \\KN \times S\\ matrix of MCMC coefficient draws of deterministic
  terms and element `lambda` contains the corresponding draws of
  inclusion parameters in case variable selection algorithms were
  employed. For time varying parameter models the coefficient matrix
  must be \\TKN \times S\\. Draws of the error covariance matrix of the
  state equation can be provided as a \\KN \times S\\ matrix in an
  additional list element.

- Sigma:

  a \\K^2 \times S\\ matrix of MCMC draws for the error
  variance-covariance matrix or a named list, where element `coeffs`
  contains a \\K^2 \times S\\ matrix of MCMC draws for the error
  variance-covariance matrix and element `lambda` contains the
  corresponding draws of inclusion parameters in case variable selection
  algorithms were employed to the covariances. For models with
  stochastic volatility the matrix must be \\TK^2 \times S\\.

- error:

  a character specifying the model that was used for the estimation of
  the covariance matrix of the error term. If `NULL` (default), it is
  inferred from the draws in argument `Sigma`, which can only
  distinguish a constant covariance matrix from a time varying one. See
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- varsel:

  a character specifying the variable selection algorithm that was
  employed. If `NULL` (default), it is inferred from the presence of
  draws of inclusion parameters. See
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- iterations:

  an integer of the number of MCMC draws in the provided coefficient
  matrices. If `NULL` (default), it is obtained from those matrices.

- burnin:

  an integer of the number of MCMC draws that were used to initialise
  the sampler. Defaults to zero, since the draws provided to the
  function are assumed to exclude them already.

## Value

An object of class `"bvarmodel"` containing the following components:

- model:

  a list containing information on the model specification.

- data:

  a list of data objects. Element `original` contains the original
  time-series objects of endogenous, unmodelled and deterministic
  variables. Element `train` contains the time-series objects `y` and
  `x` of dependent variables and regressors as well as the matrix `z` of
  regressors in SUR form.

- posterior:

  a list of posterior draws. Element `a` contains an \\S \times (K(pK +
  (1 + s)M + N) + K(K - 1) / 2)\\ "mcmc" object of the draws of the
  coefficients of the conditional mean and, if provided, the
  corresponding draws of inclusion parameters in element `lambda` and of
  the error covariance matrix of the state equation in element `sigma`.
  Element `u_sigma_inv` contains an \\S \times K^2\\ "mcmc" object of
  the draws of the inverse of the error variance-covariance matrix. For
  time varying parameter and stochastic volatility models the draws of
  all periods are appended to each other.

## Details

For the VARX model \$\$A_0 y_t = \sum\_{i = 1}^{p} A_i y\_{t-i} +
\sum\_{i = 0}^{s} B_i x\_{t - i} + C d_t + u_t\$\$ the function collects
the S draws of a Gibbs sampler (after the burn-in phase) in a
standardised object, where \\y_t\\ is a K-dimensional vector of
endogenous variables, \\A_0\\ is a \\K \times K\\ matrix of structural
coefficients. \\A_i\\ is a \\K \times K\\ coefficient matrix of lagged
endogenous variabels. \\x_t\\ is an M-dimensional vector of unmodelled,
non-deterministic variables and \\B_i\\ its corresponding coefficient
matrix. \\d_t\\ is an N-dimensional vector of deterministic terms and
\\C\\ its corresponding coefficient matrix. \\u_t\\ is an error term
with \\u_t \sim N(0, \Sigma_u)\\.

For time varying parameter and stochastic volatility models the
respective coefficients and error covariance matrix of the above model
are assumed to be time varying, respectively.

The draws of the different coefficient matrices provided in `A0`, `A`,
`B`, `C` and `Sigma` have to correspond to the same MCMC iterations.

The result is the same kind of object as the output of
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
in combination with
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md),
so it can be used with the methods of class 'bvarmodel' such as
[`irf.bvarmodel`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`predict.bvarmodel`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
or
[`summary.bvarmodel`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md).
Accordingly, the draws are stored in the parameterisation those methods
expect: the coefficients of the conditional mean are collected in a
single matrix `a` in the order in which their regressors enter the data
matrix and the draws of the error term are stored as draws of the
inverse of \\\Sigma_u\\. Draws of \\A_0\\ are reduced to the free
elements below its diagonal, which are the elements that are estimated.

Since a model that is constructed from posterior draws does not contain
the specification of its priors, the elements `priors` and `initial` of
the resulting object are empty. Use
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
and
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
if the model should be re-estimated.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r

# Get data
data("e1")
e1 <- diff(log(e1))
e1 <- window(e1, end = c(1978, 4))

# Generate model data
data <- create_bvarmodel(e1, p = 2, deterministic = "const")

# Add priors
model <- add_priors(data,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .00001))

# Set RNG seed for reproducibility
set.seed(1234567)

iterations <- 400 # Number of iterations of the Gibbs sampler
# Chosen number of iterations and burnin should be much higher.
burnin <- 100 # Number of burn-in draws
draws <- iterations + burnin # Total number of MCMC draws

y <- t(model$data$train$y)
x <- t(model$data$train$x)
tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Priors
a_mu_prior <- model$priors$a$mu # Vector of prior parameter means
a_v_i_prior <- model$priors$a$v_inv # Inverse of the prior covariance matrix

u_sigma_df_prior <- model$priors$u_sigma$df # Prior degrees of freedom
u_sigma_scale_prior <- model$priors$u_sigma$scale # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
u_sigma_i <- diag(1 / .00001, k)

# Data containers for posterior draws
draws_a <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
 # Draw conditional mean parameters
 a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)

 # Draw variance-covariance matrix
 u <- y - matrix(a, k) %*% x # Obtain residuals
 u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
 u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)

 # Store draws
 if (draw > burnin) {
  draws_a[, draw - burnin] <- a
  draws_sigma[, draw - burnin] <- solve(u_sigma_i)
 }
}

# Generate bvarmodel object
bvar_est <- bvar(data = e1,
                 y = model$data$train$y, x = model$data$train$x,
                 A = draws_a[1:18,], C = draws_a[19:21, ],
                 Sigma = draws_sigma)
```
