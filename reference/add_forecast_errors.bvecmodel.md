# Add Forecast Errors

Calculates the forecast errors of a VEC model against a test sample, in
levels.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_forecast_errors(object, test_sample, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md).

- test_sample:

  a time-series object of the endogenous variables, in levels, that
  covers the forecast periods.

- ...:

  further arguments passed to
  [`add_forecast_errors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md).

## Value

The object in `object` with `posterior$forecast_errors` added, as
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
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)
