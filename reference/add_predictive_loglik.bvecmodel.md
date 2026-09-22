# Add the Log Predictive Density of a Forecast

Scores the forecast of an object of class 'bvecmodel' against the levels
its horizon realised.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_predictive_loglik(object, test_sample = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md).

- test_sample:

  a time-series object of the endogenous variables, in levels, that
  covers the forecast periods. If `NULL` (default), the values in
  `data$test$y` of the object are used.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with `posterior$forecast$loglik` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and one column per scored period, and with the levels it was scored
against in `data$test$y`.

## Details

The forecasts of a VEC model are of the levels, so the score is of the
levels too, and it is taken in the VAR representation
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
would give: every draw is converted to its level VAR coefficients and
the realised row is evaluated against the lags of the realised rows
before it. That is the same route
[`add_forecast_errors.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md)
takes and the same one the forecast itself took, so the three describe
one model rather than three.

Everything else is as
[`add_predictive_loglik.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md)
describes it, including what each column conditions on and how a
cointegration space or a volatility that moves with time is carried
across the horizon. A model whose error correction term was scaled or
centred has to be put back with
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
first, since the score is taken against the levels themselves.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`analysis_of_stored_models`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`folder_steps`](https://franzmohr.github.io/bvartools/reference/folder_steps.md),
[`map_draws()`](https://franzmohr.github.io/bvartools/reference/map_draws.md),
[`map_models()`](https://franzmohr.github.io/bvartools/reference/map_models.md),
[`open_model()`](https://franzmohr.github.io/bvartools/reference/open_model.md),
[`open_models()`](https://franzmohr.github.io/bvartools/reference/open_models.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data("e6")
e6 <- e6 * 100
train <- window(e6, end = c(1997, 4))

# Create model
model <- create_bvecmodel(train, p = 2, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 0.0001))

model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

# Forecast the periods that were held back and score them
model <- add_forecast_input(model, n_ahead = 4)
model <- add_posterior_forecasts(model)
model <- add_predictive_loglik(model, test_sample = e6)

dim(model[["posterior"]][["forecast"]][["loglik"]])
#> [1] 20  4
```
