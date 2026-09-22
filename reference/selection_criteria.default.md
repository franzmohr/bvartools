# Selection Criteria

Calculates the criteria that a pointwise log-likelihood and a scored
forecast are enough for, whatever the model is.

## Usage

``` r
# Default S3 method
selection_criteria(object, ci = 0.95, ...)
```

## Arguments

- object:

  any object with `posterior$loglik`, `posterior$forecast$loglik` or
  both. A model of a class this package has no method for – a dynamic
  factor model of dfmtools, say – reaches this one.

- ci:

  the width of the credible bands, a value between 0 and 1. Defaults to
  0.95.

- ...:

  additional arguments.

## Value

A list of class 'selcrit', which also inherits the class of `object`,
with the element `model` and one data frame per criterion that the draws
supported, each with the columns `mean`, `median`, `qlower` and
`qupper`.
[`choose_best_model`](https://franzmohr.github.io/bvartools/reference/choose_best_model.md)
ranks a list of them and `print` shows them.

## Details

The criteria of
[`selection_criteria.bvarmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md)
fall into two groups. `AIC`, `BIC` and `HQ` charge a model for its size,
so they need a count of its free parameters, which is a property of the
model and cannot be read off its draws. `FE`, `AFE` and `RSFE` need the
forecast errors, which need the variables to be named and paired up.
Neither group is available here.

What is left needs nothing but the draws themselves:

- `LL`, the log-likelihood of the whole sample per draw, summed over the
  periods of `posterior$loglik`;

- `WAIC` and `LOOIC`, which penalise by the flexibility the fit used
  rather than by a count of parameters – see
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  for why that is what makes constant, time varying and stochastic
  volatility specifications comparable at all;

- `LPL`, the log predictive likelihood of `posterior$forecast$loglik`,
  the score of a forecast against what its horizon realised.

The periods of the `"terms"` attribute of `LPL` are numbered from one
rather than dated: this method knows nothing of where a model keeps its
sample, so it cannot say when the scored periods were. A class that can
say should write its own method, as
[`selection_criteria.bvarmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md)
does.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
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
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Any object with the draws will do, whatever produced them.
set.seed(7)
model <- list(model = list(k = 2),
              posterior = list(loglik = matrix(stats::rnorm(200, -3), 50, 4),
                               forecast = list(loglik = matrix(stats::rnorm(150, -3), 50, 3))))

criteria <- selection_criteria(model)
names(criteria)
#> [1] "model" "LL"    "WAIC"  "LOOIC" "LPL"  
criteria[["LPL"]][["mean"]]
#> [1] -7.441672
```
