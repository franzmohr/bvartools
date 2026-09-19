# Bayesian Vector Error Correction Objects

`bvec` collects the posterior draws of a vector error correction model
in an object of class 'bvecmodel'.

## Usage

``` r
bvec(
  y,
  alpha = NULL,
  beta = NULL,
  beta_x = NULL,
  beta_d = NULL,
  r = NULL,
  w = NULL,
  w_x = NULL,
  w_d = NULL,
  Gamma = NULL,
  Upsilon = NULL,
  C = NULL,
  x = NULL,
  x_x = NULL,
  x_d = NULL,
  A0 = NULL,
  Sigma = NULL,
  data = NULL,
  exogen = NULL,
  z = NULL,
  error = NULL,
  varsel = NULL,
  iterations = NULL,
  burnin = 0
)
```

## Arguments

- y:

  a time-series object of differenced endogenous variables, usually, a
  result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- alpha:

  a \\Kr \times S\\ matrix of MCMC coefficient draws of the loading
  matrix \\\alpha\\.

- beta:

  a \\Kr \times S\\ matrix of MCMC coefficient draws of cointegration
  matrix \\\beta\\ corresponding to the endogenous variables of the
  model.

- beta_x:

  a \\Mr \times S\\ matrix of MCMC coefficient draws of cointegration
  matrix \\\beta\\ corresponding to unmodelled, non-deterministic
  variables.

- beta_d:

  a \\N^{R}r \times S\\ matrix of MCMC coefficient draws of
  cointegration matrix \\\beta\\ corresponding to restricted
  deterministic terms.

- r:

  an integer of the rank of the cointegration matrix. If `NULL`
  (default), it is obtained from argument `alpha`.

- w:

  a time-series object of lagged endogenous variables in levels, which
  enter the cointegration term, usually, a result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- w_x:

  a time-series object of lagged unmodelled, non-deterministic variables
  in levels, which enter the cointegration term, usually, a result of a
  call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- w_d:

  a time-series object of deterministic terms, which enter the
  cointegration term, usually, a result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- Gamma:

  a \\(p-1)K^2 \times S\\ matrix of MCMC coefficient draws of
  differenced lagged endogenous variables or a named list, where element
  `coeffs` contains a \\(p - 1)K^2 \times S\\ matrix of MCMC coefficient
  draws of lagged differenced endogenous variables and element `lambda`
  contains the corresponding draws of inclusion parameters in case
  variable selection algorithms were employed.

- Upsilon:

  an \\sMK \times S\\ matrix of MCMC coefficient draws of differenced
  unmodelled, non-deterministic variables or a named list, where element
  `coeffs` contains a \\sMK \times S\\ matrix of MCMC coefficient draws
  of unmodelled, non-deterministic variables and element `lambda`
  contains the corresponding draws of inclusion parameters in case
  variable selection algorithms were employed.

- C:

  an \\KN^{UR} \times S\\ matrix of MCMC coefficient draws of
  unrestricted deterministic terms or a named list, where element
  `coeffs` contains a \\KN^{UR} \times S\\ matrix of MCMC coefficient
  draws of deterministic terms and element `lambda` contains the
  corresponding draws of inclusion parameters in case variable selection
  algorithms were employed.

- x:

  a time-series object of \\K(p - 1)\\ differenced endogenous variables.

- x_x:

  a time-series object of \\Ms\\ differenced unmodelled regressors.

- x_d:

  a time-series object of \\N^{UR}\\ deterministic terms that do not
  enter the cointegration term.

- A0:

  either a \\K^2 \times S\\ matrix of MCMC coefficient draws of
  structural parameters or a named list, where element `coeffs` contains
  a \\K^2 \times S\\ matrix of MCMC coefficient draws of structural
  parameters and element `lambda` contains the corresponding draws of
  inclusion parameters in case variable selection algorithms were
  employed.

- Sigma:

  a \\K^2 \times S\\ matrix of MCMC draws for the error
  variance-covariance matrix or a named list, where element `coeffs`
  contains a \\K^2 \times S\\ matrix of MCMC draws for the error
  variance-covariance matrix and element `lambda` contains the
  corresponding draws of inclusion parameters in case variable selection
  algorithms were employed to the covariances. For models with
  stochastic volatility the matrix must be \\TK^2 \times S\\.

- data:

  the original time-series object of endogenous variables in levels. If
  `NULL` (default), it is reconstructed from arguments `y` and `w`.

- exogen:

  the original time-series object of unmodelled variables in levels. If
  `NULL` (default), it is reconstructed from arguments `w_x` and `x_x`.

- z:

  a \\TK \times K (r + (p-1)K + sM + N)\\ data matrix, with \\K(K -
  1)/2\\ further columns for the contemporaneous endogenous variables of
  a structural model, usually, a result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).
  If `NULL` (default), it is generated from the provided data.

- error:

  a character specifying the model that was used for the estimation of
  the covariance matrix of the error term. If `NULL` (default), it is
  inferred from the draws in argument `Sigma`, which can only
  distinguish a constant covariance matrix from a time varying one. See
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- varsel:

  a character specifying the variable selection algorithm that was
  employed. If `NULL` (default), it is inferred from the presence of
  draws of inclusion parameters. See
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- iterations:

  an integer of the number of MCMC draws in the provided coefficient
  matrices. If `NULL` (default), it is obtained from those matrices.

- burnin:

  an integer of the number of MCMC draws that were used to initialise
  the sampler. Defaults to zero, since the draws provided to the
  function are assumed to exclude them already.

## Value

An object of class `"bvecmodel"` containing the following components:

- model:

  a list containing information on the model specification.

- data:

  a list of data objects. Element `original` contains the original
  time-series objects of endogenous and unmodelled variables in levels.
  Element `train` contains the time-series objects `y` of differenced
  endogenous variables, `w` of the regressors of the error correction
  term and `x` of the remaining regressors as well as the matrix `z` of
  regressors in SUR form.

- posterior:

  a list of posterior draws. Element `a` contains an \\S \times (Kr +
  K((p - 1)K + sM + N^{UR}) + K(K - 1) / 2)\\ "mcmc" object of the draws
  of the loading matrix and the coefficients of the remaining regressors
  and, if provided, the corresponding draws of inclusion parameters in
  element `lambda` and of the error covariance matrix of the state
  equation in element `sigma`. Element `beta` contains an \\S \times
  ((K + M + N^{R})r)\\ "mcmc" object of the draws of the cointegration
  matrix. Element `u_sigma_inv` contains an \\S \times K^2\\ "mcmc"
  object of the draws of the inverse of the error variance-covariance
  matrix.

## Details

For the vector error correction model with unmodelled exogenous
variables (VECX) \$\$A_0 \Delta y_t = \alpha \beta^\prime
\begin{pmatrix} y\_{t-1} \\ x\_{t-1} \\ d^{R}\_{t-1} \end{pmatrix} +
\sum\_{i = 1}^{p-1} \Gamma_i \Delta y\_{t-i} + \sum\_{i = 0}^{s-1}
\Upsilon_i \Delta x\_{t-i} + C^{UR} d^{UR}\_t + u_t\$\$ the function
collects the \\S\\ draws of a Gibbs sampler in a standardised object,
where \\\Delta y_t\\ is a K-dimensional vector of differenced endogenous
variables and \\A_0\\ is a \\K \times K\\ matrix of structural
coefficients. \\\alpha\\ is the \\K \times r\\ loading matrix and
\\\beta\\ the \\(K + M + N^{R}) \times r\\ cointegration matrix of the
error correction term, where \\y\_{t-1}\\, \\x\_{t-1}\\ and
\\d^{R}\_{t-1}\\ are the first lags of endogenous, exogenous variables
in levels and restricted deterministic terms, respectively. \\\Gamma_i\\
is a coefficient matrix of lagged differenced endogenous variabels.
\\\Delta x_t\\ is an M-dimensional vector of unmodelled,
non-deterministic variables and \\\Upsilon_i\\ its corresponding
coefficient matrix. \\d_t\\ is an \\N^{UR}\\-dimensional vector of
unrestricted deterministics and \\C^{UR}\\ the corresponding coefficient
matrix. \\u_t\\ is an error term with \\u_t \sim N(0, \Sigma_u)\\.

For time varying parameter and stochastic volatility models the
respective coefficients and error covariance matrix of the above model
are assumed to be time varying, respectively.

The draws of the different coefficient matrices provided in `alpha`,
`beta`, `beta_x`, `beta_d`, `A0`, `Gamma`, `Upsilon`, `C` and `Sigma`
have to correspond to the same MCMC iteration.

The result is the same kind of object as the output of
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
in combination with
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md),
so it can be used with the methods of class 'bvecmodel' such as
[`vec_to_var.bvecmodel`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)
or
[`summary.bvecmodel`](https://franzmohr.github.io/bvartools/reference/summary.bvecmodel.md).
Accordingly, the draws are stored in the parameterisation those methods
expect. The error correction term is described by \\\alpha\\ and
\\\beta\\ and not by \\\Pi = \alpha \beta^\prime\\, since the latter
cannot be decomposed into the former. The draws of \\\alpha\\ are
collected in the same matrix as the coefficients of the remaining
regressors and the draws of \\\beta\\ of all cointegration relations are
collected in a matrix of their own. Draws of the error term are stored
as draws of the inverse of \\\Sigma_u\\ and draws of \\A_0\\ are reduced
to the free elements below its diagonal, which are the elements that are
estimated.

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
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r

# Load data
data("e6")
# Generate model
model <- create_bvecmodel(e6, p = 4, r = 1, const = "unrestricted", seasonal = "unrestricted")
# Obtain data matrices
y <- t(model$data$train$y)
w <- t(model$data$train$w)
x <- t(model$data$train$x)

# Reset random number generator for reproducibility
set.seed(1234567)

iterations <- 400 # Number of iterations of the Gibbs sampler
# Chosen number of iterations should be much higher, e.g. 30000.

burnin <- 100 # Number of burn-in draws
draws <- iterations + burnin

r <- 1 # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta
k_gamma <- k * k_x

# Set uninformative priors
a_mu_prior <- matrix(0, k_alpha + k_gamma) # Vector of prior parameter means
a_v_i_prior <- diag(0, k_alpha + k_gamma) # Inverse of the prior covariance matrix

v_i <- 0
p_tau_i <- diag(1, k_w)

u_sigma_df_prior <- r # Prior degrees of freedom
u_sigma_scale_prior <- diag(0.01, k) # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(c(1, -4), k_w, r)
u_sigma_i <- diag(1 / .0001, k)
g_i <- u_sigma_i

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, sigma_i = u_sigma_i,
                         v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                         x = x,
                         gamma_mu_prior = a_mu_prior,
                         gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma

  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  u_sigma_scale_post <- solve(tcrossprod(u) +
     v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  u_sigma <- solve(u_sigma_i)

  # Update g_i
  g_i <- u_sigma_i

  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- u_sigma
  }
}

# Number of non-deterministic coefficients
k_nondet <- (k_x - 4) * k

# Generate bvecmodel object
bvec_est <- bvec(y = model$data$train$y, w = model$data$train$w,
                 x = model$data$train$x[, 1:6],
                 x_d = model$data$train$x[, 7:10],
                 r = r,
                 alpha = draws_alpha,
                 beta = draws_beta,
                 Gamma = draws_gamma[1:k_nondet,],
                 C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
                 Sigma = draws_sigma)
```
