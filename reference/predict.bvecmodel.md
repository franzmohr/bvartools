# Predict Method for Objects of Class bvecmodel

Summarises the forecasts of a VEC model, in levels.

## Usage

``` r
# S3 method for class 'bvecmodel'
predict(object, n_ahead = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html).

- n_ahead:

  number of steps ahead at which to predict. If `NULL` (default), every
  period simulated by
  [`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html).

- ...:

  additional arguments passed to
  [`predict.bvarmodel`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md).

## Value

A time-series object of class `"bvarprd"`, as returned by
[`predict.bvarmodel`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md).

## Details

The forecasts are of the levels, and so is the history they are shown
with: it is recovered from the differences and the error correction term
the way
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
recovers it.

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
[`add_seed.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_seed.bvarmodel.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md)

## Examples

``` r

data("e6")
model <- create_bvecmodel(e6 * 100, p = 1, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 4))
prd <- predict(model)
plot(prd)


```
