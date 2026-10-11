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
  [`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html).

- n_ahead:

  an integer of the forecast horizon.

- deterministic, exogen:

  the values of the deterministic terms and of the unmodelled variables
  in the forecast periods, as for
  [`add_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md).

- ...:

  further arguments passed to
  [`prepare_forecast_input`](https://rdrr.io/pkg/bayests/man/prepare_forecast_input.html).

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

The forecast starts right after the estimation sample. Argument
`y_start` of
[`add_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md)
is refused here; to start a forecast from values after the sample,
convert the model with
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
and give them, in levels, to the resulting 'bvarmodel'.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_seed.bvarmodel.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

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

model <- add_forecast_input(model, n_ahead = 4)
model <- add_posterior_forecasts(model)
plot(predict(model))


```
