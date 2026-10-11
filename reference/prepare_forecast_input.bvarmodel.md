# Prepare Forecast Input

Generates data matrices serving as input for forecasting simulation for
objects of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
prepare_forecast_input(
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
  of the model. Required if the model has such variables. See 'Details'.

- y_start:

  an optional time-series object with the values of the endogenous
  variables the forecast starts from. The forecast periods are the
  `n_ahead` periods after its last observation, which may lie after the
  end of the estimation sample. If `NULL` (default), the forecast starts
  after the end of the estimation sample. See 'Details'.

- ...:

  additional arguments.

## Value

A list with elements `h`, the forecast horizon, and `x`, the
out-of-sample regressors: `h` rows, one per period, by one column per
regressor. That is the compact layout, the same one a coefficient matrix
is `k` by; the SUR layout this used to return spread every regressor
over `k` columns and was `k^2` the size for no extra content. The lags
of the endogenous variables that the estimation sample does not reach
are the forecasts of earlier periods, which the forecast fills in as it
goes; they are zero here, since the samplers refuse a missing value
anywhere in `x`.

If `y_start` is given, the list has two more elements: `start`, the time
of the first forecast period, and `y`, the observations of `y_start`
after the end of the estimation sample together with the last `p` the
forecast starts from, a time-series object with one column per
endogenous variable.

## Details

The regressors of a forecast period contain the values of the
unmodelled, non-deterministic variables in that period and in the `s`
periods before it. Argument `exogen` therefore has to cover the last `s`
periods of the estimation sample as well as the `n_ahead` forecast
periods, at the frequency of the model. A series that starts with the
first forecast period, or ends before the last one, is refused with a
message that names the periods it has to cover.

If `deterministic` is not given, the deterministic terms are continued
from the estimation sample, which works for a constant, a linear trend,
seasonal dummies and the dummy variables of
[`add_dummy_variables`](https://franzmohr.github.io/bvartools/reference/add_dummy_variables.md).

**Starting a forecast after the end of the estimation sample.** By
default the first forecast period is the one after the last period of
the estimation sample, and its lagged regressors are the last `p`
observations of that sample. Argument `y_start` moves the start: the
first forecast period is then the one after the last observation of
`y_start`, and the lags are its last `p` observations. This is how a
forecast starts from a value that was not part of the estimation, such
as a nowcast of the current month or a first release that came in after
the model was estimated, without estimating the model again. `y_start`
holds the endogenous variables of the model, by name and in the units of
the estimation sample, at the frequency of the model. It has to end in
the last period of the estimation sample or after it, and its last `p`
observations must not be missing; values before those are not read. When
it ends in the last period of the estimation sample and holds the same
values, the result is the same as without it.

The other inputs follow the new forecast periods: `exogen` has to cover
the `s` periods before the first of them and the `n_ahead` periods
themselves, and `deterministic` the forecast periods. The automatic
continuation of the deterministic terms counts through the periods
between the end of the sample and the start of the forecast, so that a
trend keeps counting and seasonal dummies fall on the calendar months or
quarters of the forecast periods.

The coefficients, and for models with time varying coefficients or
stochastic volatility the states, are those estimated on the sample.
Nothing is learnt from the start values beyond their role as lags. Some
uses of a forecast therefore need it to start right after the sample:
see
[`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html)
and
[`add_predictive_loglik`](https://rdrr.io/pkg/bayests/man/add_predictive_loglik.html).
A model estimated with `missing = "estimate"` starts each draw's
forecast from its own completed panel and does not take `y_start`.

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

# Generate forcast input
fcst_input <- prepare_forecast_input(model, n_ahead = 4)

# Start the forecast two quarters after the end of the sample, from
# values the model was not estimated with
short <- create_bvarmodel(window(e1, end = c(1978, 4)), p = 2,
                          deterministic = "const",
                          iterations = 20, burnin = 10)
fcst_input <- prepare_forecast_input(short, n_ahead = 4,
                                     y_start = window(e1, end = c(1979, 2)))
fcst_input[["start"]]
#> [1] 1979.5

```
