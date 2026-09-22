# Seed of the Posterior Simulation

Sets the seed with which the posterior of a model is simulated.

## Usage

``` r
add_seed(object, seed, ...)

# S3 method for class 'bvarmodel'
add_seed(object, seed, ...)

# S3 method for class 'bvecmodel'
add_seed(object, seed, ...)

# S3 method for class 'modellist'
add_seed(object, seed, ...)

# S3 method for class 'expandingwindow'
add_seed(object, seed, ...)
```

## Arguments

- object:

  a model, at any point before the posterior is simulated: an object of
  class 'bvarmodel' or 'bvecmodel', a model of another package that
  provides an `add_seed` method for it, or a list of such models of
  class 'modellist' or 'expandingwindow'.

- seed:

  a non-negative whole number no larger than `.Machine$integer.max`.

- ...:

  further arguments passed to or from other methods.

## Value

`object` with the seed set.

## Details

Calling this function is optional.
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
already stores a seed, drawn from R's random number generator, in
element `seed` of `object$model`, unless the model has one. `add_seed`
sets that element, for example to give a model a seed that depends
neither on the state of R's generator nor on the worker of a cluster
that happens to simulate it. It can be called before
`add_initial_values`, which then keeps the seed it is given, or after
it, which replaces the drawn one.

The seed is part of the model.
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
writes it as attribute `seed` of group `/model`, where the BayesTS
executable reads it, and
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
reads it back. The internal samplers of
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
draw with it too, so a model with a given seed gives the same draws
however R's generator stands. BayesTS and this package use different
generators, so the same seed gives different draws in the two.

A list of models gets the seeds `seed`, `seed + 1`, ..., one per model
in the order of its elements, counting through nested lists, so that no
two of its models draw the same numbers. A model is any element with an
`add_seed` method of its own, so the models of other packages that
provide one, such as the dynamic factor models of dfmtools, are counted
with the rest. Elements that are not estimated, such as external
forecasts, and elements without a method are left as they are and are
not counted.

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

# Generate model
model <- create_bvarmodel(data = e1, p = 2, iterations = 100, burnin = 50)

# Set the seed of the posterior simulation
model <- add_seed(model, 20260916)
model[["model"]][["seed"]]
#> [1] 20260916
```
