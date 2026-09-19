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
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md)

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
#>    invest 1   4.103   3.699  5.1893
#>    income 1   1.618   1.395  1.3465
#>      cons 1   1.438   1.494  0.9281
#>    invest 2   3.447   6.021  4.4958
#>    income 2   1.337   1.448  1.1919
#>      cons 2   1.230   1.038  0.7287
#> 
#> 
#> Root mean squared forecast errors (RMSFE)
#> 
#>  Variable h Model 1 Model 2 Model 3
#>    invest 1   4.988   4.512  6.1921
#>    income 1   1.901   1.760  1.6416
#>      cons 1   1.701   1.722  1.2042
#>    invest 2   4.121   7.015  5.3257
#>    income 2   1.892   1.843  1.4546
#>      cons 2   1.536   1.405  0.9832
```
