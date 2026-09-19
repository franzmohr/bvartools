# Align Observations Across Models

Restricts each model in an object of class 'modellist' to the set of
observations common to all models, ensuring that comparisons are
computed on the same underlying sample.

## Usage

``` r
# S3 method for class 'modellist'
align_model_obs(object, ...)
```

## Arguments

- object:

  an object of class 'modellist'.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'modellist'.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data(e6)

# Create multiple VAR models
model_1 <- create_bvarmodel(diff(e6), p = 0:2, deterministic = "const")

# Create multiple VEC models
model_2 <- create_bvecmodel(e6, p = 1, r = 0:1, const = "unrestricted")

# Combine the models in one object
model <- combine_models(model_1, model_2)

model <- align_model_obs(model)
```
