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

- ...:

  additional arguments.

## Value

A list with elements `h`, the forecast horizon, and `x`, the
out-of-sample regressors: `h` rows, one per period, by one column per
regressor. That is the compact layout, the same one a coefficient matrix
is `k` by; the SUR layout this used to return spread every regressor
over `k` columns and was `k^2` the size for no extra content.

## Details

The regressors of a forecast period contain the values of the
unmodelled, non-deterministic variables in that period and in the `s`
periods before it. Argument `exogen` therefore has to cover the last `s`
periods of the estimation sample as well as the `n_ahead` forecast
periods, at the frequency of the model. A series that starts with the
first forecast period, or ends before the last one, is refused with a
message that names the periods it has to cover.

If `deterministic` is not given, the deterministic terms are continued
from the estimation sample, which works for a constant, a linear trend
and seasonal dummies.

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

```
