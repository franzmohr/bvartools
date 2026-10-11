# Add Forecasts

Calculates and adds forecasts to an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_posterior_forecasts(
  object,
  forecast_states = NULL,
  scenario = NULL,
  forecast_quantile = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvarmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
  and
  [`add_forecast_input`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html).

- forecast_states:

  character, what a model with time varying coefficients or stochastic
  volatility does with them over the forecast horizon. `"simulate"`
  carries each draw's random walks forward, one step per period, so that
  the forecasts are draws from the posterior predictive distribution of
  the estimated model. `"hold"` keeps the coefficients and volatilities
  at their values in the last sample period, which gives forecasts
  conditional on no further drift and narrower intervals, and is what
  earlier versions of the package did. If `NULL` (default), the value in
  `object$model$forecast_states` is used, and `"simulate"` when there is
  none. Models with constant coefficients and volatility are unaffected.

- scenario:

  an optional data frame with columns `period`, `variable` and `value`,
  one row per value a variable is held at in a forecast period: a
  conditional forecast. `period` counts forecast periods from one and
  `variable` is a name or a position among the endogenous variables. It
  is stored in `data$forecast$constraints`, which is where the forecast
  reads it from; set that element to `NULL` to forecast without it
  again. Available for models with constant coefficients and
  `error = "wishart"` or `"gamma"`, and for a grid of quantiles.

- forecast_quantile:

  for a grid of quantiles (see `quantile_grid` in
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)),
  a level in \\(0, 1)\\ every draw of the forecast is taken at, which
  gives quantile paths rather than draws from the predictive
  distribution: the paths of Chavleishvili and Manganelli (2019), whose
  difference with and without a `scenario` is an impulse response at
  that quantile. If `NULL` (default), the value in
  `object$model$forecast_quantile` is used, and the levels are drawn at
  random when there is none. Zero removes a stored one.

- ...:

  arguments passed forward to method.

## Value

The object in `object` with `posterior$forecast$forecasts` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and \\Kh\\ columns, stacked by period: the \\K\\ variables of the
first forecast period, then those of the second, and so on.
`posterior$forecast` is the group everything the forecast periods
produce hangs below, the `errors` of
[`add_forecast_errors`](https://rdrr.io/pkg/bayests/man/add_forecast_errors.html)
beside these.
[`predict`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
summarises them. A `forecast_states` that was given is stored in
`model$forecast_states`.

Simulating the volatility forward needs the variance of the
log-volatility innovations, `posterior$u_sigma_inv$sigma`, which
[`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
stores. A stochastic volatility model fitted with an earlier version of
the package lacks it and stops with an error unless
`forecast_states = "hold"`.

## Details

A model with a seed, see
[`add_seed`](https://rdrr.io/pkg/bayests/man/add_seed.html), draws its
forecast from a seed derived from its own, so two forecasts of it give
the same draws however R's generator stands. A model without one draws
from R's generator as it stands.

A forecast that starts after the end of the estimation sample, from
argument `y_start` of
[`add_forecast_input`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html),
takes its lags from the start values and its coefficients from the
sample. **Three things that start from the last period of the sample are
refused for it**: a `scenario`, which is drawn from the end of the
sample; `forecast_states = "simulate"` for a model with time varying
coefficients or stochastic volatility, which would not step the states
through the periods up to the start values (`"hold"` is available); and
the discounted model, for the same reason.

## References

Chavleishvili, S., & Manganelli, S. (2019). Forecasting and stress
testing with quantile vector autoregression. *ECB Working Paper*, 2330.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
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

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Add data used for forecast calculation
model <- add_forecast_input(model, n_ahead = 4)

# Add forecasts
model <- add_posterior_forecasts(model)

```
