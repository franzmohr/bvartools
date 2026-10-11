# Seed of the Posterior Simulation of a VAR or VEC Model

The methods of
[`add_seed`](https://rdrr.io/pkg/bayests/man/add_seed.html) for VAR and
VEC models: they set the seed with which the posterior is simulated, in
element `seed` of `object$model`.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_seed(object, seed, ...)

# S3 method for class 'bvecmodel'
add_seed(object, seed, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel' or 'bvecmodel'.

- seed:

  a non-negative whole number no larger than `.Machine$integer.max`.

- ...:

  not used.

## Value

`object` with the seed set.

## Details

[`write_to_hdf5`](https://rdrr.io/pkg/bayests/man/write_to_hdf5.html)
writes the seed as attribute `seed` of group `/model`, where the BayesTS
executable reads it, and the internal samplers of
[`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
draw with it too, so a model with a given seed gives the same draws
however R's generator stands. So do the simulations that follow the
coefficients:
[`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html),
scenarios included,
[`add_predictive_loglik`](https://rdrr.io/pkg/bayests/man/add_predictive_loglik.html),
[`add_sign_restrictions`](https://rdrr.io/pkg/bayests/man/add_sign_restrictions.html)
and
[`add_sign_zero_restrictions`](https://rdrr.io/pkg/bayests/man/add_sign_zero_restrictions.html).
Each of these three stages – the forecast, the predictive score and the
rotations – draws from a seed of its own derived from the model's, as
the BayesTS command line seeds its stages from `/model/seed`, so that a
forecast repeated gives the same draws and does not reuse the random
numbers the chain started from. A model without a seed draws from R's
generator as it stands throughout. BayesTS and this package use
different generators, so the same seed gives different draws in the two.
A list of models is numbered by the method of bayests.

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
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r
# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model
model <- create_bvarmodel(data = e1, p = 2, iterations = 100, burnin = 50)

# Set the seed of the posterior simulation
model <- add_seed(model, 20260916)
model[["model"]][["seed"]]
#> [1] 20260916
```
