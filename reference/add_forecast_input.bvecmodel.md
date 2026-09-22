# Add Forecast Input Data

Prepares the input data of the forecast periods of a VEC model.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_forecast_input(
  object,
  n_ahead = 10,
  deterministic = NULL,
  exogen = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- n_ahead:

  an integer of the forecast horizon.

- deterministic, exogen:

  the values of the deterministic terms and of the unmodelled variables
  in the forecast periods, as for
  [`add_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md).

- ...:

  further arguments passed to
  [`prepare_forecast_input`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.md).

## Value

The object in `object` with the forecast horizon in `model$h` and the
regressors of the forecast periods, in levels, in `data$forecast$x`.

## Details

A VEC model is forecast in levels: it is the same model as its VAR
representation, and the forecast of the differences follows from that of
the levels. The regressors of the forecast periods are therefore those
of the VAR in levels, and are assembled from the data and specification
that
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
builds. The posterior draws are not converted.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
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
