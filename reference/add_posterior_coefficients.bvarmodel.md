# Posterior Simulation of Model Coefficients

Forwards model input to posterior simulation functions for vector
autoregressive models.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_posterior_coefficients(
  object,
  posterior_function = NULL,
  chains = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvarmodel', usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  in combination with
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
  and
  [`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md).

- posterior_function:

  the function to be applied to the model in argument `object`. If
  `NULL` (default), internal functions are used.

- chains:

  the number of chains to simulate. If `NULL` (default), the value in
  `object$model$chains` is used, and one chain when there is none. Each
  chain is the same simulation with a seed of its own – the first with
  the seed of the model, so that one chain draws what it always has –
  and the chains are pooled, one after the other, in the draws of
  `posterior`, so that every later step uses all of them. Their number,
  when above one, is stored in `object$model$chains`, and
  [`chain_diagnostics`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md)
  compares them. A `posterior_function` is called once per chain. Not
  available for discounted models, whose posterior is not a chain.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `posterior` added. Each of its
elements is a list whose element `coeffs` holds the draws after burn-in
as a [`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one
row per draw and one column per parameter:

- `a`:

  the coefficients, \\M\\ columns or \\TM\\ for TVP models, in the order
  of the columns of `data$train$z` with the contemporaneous coefficients
  of structural models last. With variable selection element `lambda`
  holds the inclusion indicators, and for TVP models element `sigma` the
  state variances.

- `u_sigma_inv`:

  the inverse error covariance matrix, \\K^2\\ columns, or \\TK^2\\ if
  the error variances vary over time.

- `u_omega_inv`:

  for all errors but `"wishart"`, the error precisions, \\K\\ columns or
  \\TK\\.

- `psi`:

  for `"gamma+covar"` and `"sv+covar"`, the error covariance
  coefficients.

- `u_scale`:

  for `error = "ald"`, the \\K\\ scales of the asymmetric Laplace
  distribution.

Elements that do not apply to a model are absent or `NULL`. Note that
[`bvar`](https://franzmohr.github.io/bvartools/reference/bvar.md)
expects draws in the transposed orientation.

## Details

Unless `posterior_function` is specified, the function forwards the
model input to the package's own posterior functions.

The internal samplers draw with the seed in `object$model$seed`, which
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
sets and
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
replaces. R's random number generator is set to that seed, with R's
default kinds, for the simulation and put back as it was afterwards. A
call of [`set.seed()`](https://rdrr.io/r/base/Random.html) between
[`add_initial_values()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
and this function therefore does not change the draws. A model without a
seed draws from R's generator as it stands. A `posterior_function` is
called as it is and decides itself what to do with the seed; see
[`bayests_posterior`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md).

A sampler that cannot run raises its error rather than returning
something. The message names what about the input it could not work
with. Applied to a list of models – a 'modellist', an 'expandingwindow'
or, in bgvars, a 'gvarmodel' – that ends the run on the first
specification that fails, rather than leaving that one without a
posterior and carrying it into whatever reads the results.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_files()`](https://franzmohr.github.io/bvartools/reference/bayests_files.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`chain_diagnostics()`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)
```
