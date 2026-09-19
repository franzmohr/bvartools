# Add Forecasts

Simulates forecasts of a VEC model, in levels, from its posterior draws.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_posterior_forecasts(object, forecast_states = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  and
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md).

- forecast_states:

  character, what a VEC model with time varying coefficients or
  stochastic volatility does with them over the forecast horizon.
  `"simulate"` carries each draw's random walks forward, one step per
  period – the loadings and short-run coefficients, the cointegration
  vectors through their state equation, the covariance block and the
  log-volatilities – and rebuilds the VAR in levels from them at every
  period, so that the forecasts are draws from the posterior predictive
  distribution of the estimated model. `"hold"` keeps them at their
  values in the last sample period. If `NULL` (default), the value in
  `object$model$forecast_states` is used, and `"simulate"` when there is
  none. Models with constant coefficients and volatility are unaffected.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with `posterior$forecast` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and \\Kh\\ columns of forecasts of the levels, stacked by period.
`predict` summarises them. A `forecast_states` that was given is stored
in `model$forecast_states`.

## Details

The forecasts are those of the VAR in levels the VEC model implies. For
a model with constant coefficients they are the forecasts
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
followed by
[`add_posterior_forecasts.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md)
gives. For a model with time varying coefficients they are not: its VAR
representation has no state equation of its own and holds its
coefficients at the last period, while this method steps the states of
the VEC model and converts them to levels anew in every forecast period.

Simulating the volatility forward needs `posterior$u_sigma_inv$sigma`,
the variance of the log-volatility innovations, which
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
stores. A model with stochastic volatility estimated with an earlier
version of the package lacks it and stops with an error unless
`forecast_states = "hold"`.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)
