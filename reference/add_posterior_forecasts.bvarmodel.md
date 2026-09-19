# Add Forecasts

Calculates and adds forecasts to an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_posterior_forecasts(object, forecast_states = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  and
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md).

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

- ...:

  arguments passed forward to method.

## Value

The object in `object` with `posterior$forecast` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and \\Kh\\ columns, stacked by period: the \\K\\ variables of the
first forecast period, then those of the second, and so on.
[`predict`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
summarises them. A `forecast_states` that was given is stored in
`model$forecast_states`.

Simulating the volatility forward needs the variance of the
log-volatility innovations, `posterior$u_sigma_inv$sigma`, which
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
stores. A stochastic volatility model fitted with an earlier version of
the package lacks it and stops with an error unless
`forecast_states = "hold"`.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
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
