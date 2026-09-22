# Model Selection Criteria

Calculates model selection criteria for a list of Bayesian models.

## Usage

``` r
# S3 method for class 'modellist'
selection_criteria(object, ...)
```

## Arguments

- object:

  an object of class 'modellist'.

- ...:

  further arguments passed to or from other methods.

## Value

A list of class 'selcritlist' with one object of class 'selcrit' per
model, as described in
[`selection_criteria.bvarmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
which
[`choose_best_model`](https://franzmohr.github.io/bvartools/reference/choose_best_model.md)
compares.

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
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md)

## Examples

``` r

# Load data
data("e1")
orig <- diff(log(e1)) * 100
train <- window(orig, end = c(1982, 2))


# Create model
model <- create_bvarmodel(data = train, p = 0:2, deterministic = "const",
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

# Add forecast errors
model <- add_forecast_errors(model, test_sample = orig)

# Calculate selection criteria
sel <- selection_criteria(model)
sel
#> 
#> 
#> ------------------------------------------
#> Out-of-sample
#> ------------------------------------------
#> 
#> Mean absolute forecast errors (MAFE)
#> 
#>  Variable h Model 1 Model 2 Model 3
#>    invest 1   2.288   4.254   4.069
#>    income 1   1.868   1.648   1.629
#>      cons 1   1.619   1.141   1.072
#>    invest 2   3.986   4.959   4.911
#>    income 2   1.346   1.173   1.471
#>      cons 2   1.147   1.309   1.013
#> 
#> 
#> Root mean squared forecast errors (RMSFE)
#> 
#>  Variable h Model 1 Model 2 Model 3
#>    invest 1   2.836   5.074   4.871
#>    income 1   2.156   2.019   1.952
#>      cons 1   1.889   1.404   1.346
#>    invest 2   4.713   5.727   6.295
#>    income 2   1.593   1.390   1.742
#>      cons 2   1.307   1.620   1.169
```
