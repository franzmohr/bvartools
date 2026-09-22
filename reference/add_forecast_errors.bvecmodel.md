# Add Forecast Errors

Calculates the forecast errors of a VEC model against a test sample, in
levels.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_forecast_errors(object, test_sample = NULL, ...)
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

  further arguments passed to
  [`add_forecast_errors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md).

## Value

The object in `object` with `posterior$forecast$errors` added, as
described in
[`add_forecast_errors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md).

## Details

The forecasts of a VEC model are of the levels, so the errors are taken
against the levels of `test_sample`, exactly as for the VAR
representation
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
would give.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
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
