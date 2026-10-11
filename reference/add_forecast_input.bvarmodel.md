# Add Forecast Input Data

Generates and adds data matrices for forecast simulation to the elements
of an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_forecast_input(
  object,
  n_ahead = 10,
  deterministic = NULL,
  exogen = NULL,
  y_start = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvarmodel'.

- n_ahead:

  number of steps ahead at which to predict.

- deterministic:

  a time-series object with deterministic data. If not specified, the
  function will try to identify the deterministic terms automatically.
  If this is not successful, an error message we be returned.

- exogen:

  a time-series object with the unmodelled, non-deterministic variables
  of the model. Required if the model has such variables. It has to
  cover the `s` periods before the first forecast period, which are the
  last `s` periods of the estimation sample, as well as the `n_ahead`
  forecast periods, because the regressors of a forecast period include
  the lags of these variables. See
  [`prepare_forecast_input`](https://rdrr.io/pkg/bayests/man/prepare_forecast_input.html).

- y_start:

  an optional time-series object with the endogenous variables of the
  model, by name and in the units of the estimation sample, that the
  forecast starts from: the forecast periods are the `n_ahead` periods
  after its last observation, and the lags of the first one are its last
  `p` observations. It may end after the estimation sample, so that a
  forecast can start from a nowcast or from data that came in after the
  model was estimated; `exogen` and `deterministic` then have to cover
  the new forecast periods. If `NULL` (default), the forecast starts
  after the end of the estimation sample. See 'Details' in
  [`prepare_forecast_input`](https://rdrr.io/pkg/bayests/man/prepare_forecast_input.html).

- ...:

  arguments passed forward to method.

## Value

The object in `object` with `model$h` set to `n_ahead` and
`data$forecast$x` added, the regressors of the forecast periods in the
layout of `data$train$x`, one row per period. If `y_start` is given,
`model$forecast_start` is the time of the first forecast period, which
[`predict`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
dates the forecasts by, and `data$forecast$y` holds the start values
after the end of the sample, which it shows as the history the forecasts
continue.

## See also

Other posterior simulation:
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
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

# Add data used for forecast calculation
model <- add_forecast_input(model, n_ahead = 4)

# A model estimated up to 1978Q4 whose forecast starts from the values
# of 1979Q1, which it was not estimated with
short <- create_bvarmodel(window(e1, end = c(1978, 4)), p = 2,
                          deterministic = "const",
                          iterations = 20, burnin = 10)
short <- add_priors(short,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))
short <- add_forecast_input(short, n_ahead = 4,
                            y_start = window(e1, end = c(1979, 1)))
short[["model"]][["forecast_start"]]
#> [1] 1979.25
```
