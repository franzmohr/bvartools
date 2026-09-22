# Model Selection Criteria

Calculates model selection criteria for an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
selection_criteria(object, ci = 0.95, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- ...:

  further arguments passed to or from other methods.

## Value

A list of class 'selcrit', which also inherits the class of the model,
with the element `model` and one data frame per criterion. If the model
contains `posterior$loglik`, these are `LL`, `AIC`, `BIC`, `HQ`, `WAIC`
and `LOOIC`, each with the columns `mean`, `median`, `qlower` and
`qupper`, where bands that do not apply are `NA`. If it contains
`posterior$forecast$errors`, these are `FE`, `AFE` and `RSFE`, the
forecast errors and their absolute and root squared values, with the
columns `variable`, `h`, `mean`, `median`, `qlower` and `qupper`. For
`RSFE` only `mean` is a statistic of its own, the root of the mean
squared error. Its median and quantiles are those of the squared errors,
rooted, and so differ from `AFE`'s only by the interpolation between
neighbouring draws that
[`quantile`](https://rdrr.io/r/stats/quantile.html) and
[`median`](https://rdrr.io/r/stats/median.html) do – which, the root of
an average of squares being at least the average of the roots, can only
raise them. They add nothing to `AFE`'s. If it contains
`posterior$forecast$loglik`, the score of its forecast against
`data$test$y`, there is `LPL`, the log predictive likelihood: the log of
the mean of the draws of each scored period's predictive density, summed
over the periods. Attribute `"terms"` holds those per period and `"nse"`
their numerical standard error.

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
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 1:2, deterministic = "const",
                          iterations = 10, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Add log-likelihoods
model <- add_posterior_loglik(model)

# Calculate selection criteria
sel <- selection_criteria(model)
sel
#> 
#> ------------------------------------------
#> In-sample
#> ------------------------------------------
#> 
#> Log-likelihood
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -511.2 -511.8          -515.5           -507.8
#>  Model 2 -499.1 -499.2          -501.8           -495.5
#> 
#> 
#> Akaike Information Criterion (AIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1040   1040                                 
#>  Model 2 1031   1031                                 
#> 
#> 
#> Bayesian Information Criterion (BIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1085   1085                                 
#>  Model 2 1098   1098                                 
#> 
#> 
#> Hannan-Quinn Criterion (HQ)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1058   1058                                 
#>  Model 2 1058   1058                                 
#> 
#> 
#> Widely Applicable Information Criterion (WAIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1047   1047           980.3             1114
#>  Model 2 1024   1024           959.7             1089
#> 
#> Periods with a pointwise log-likelihood variance above 0.4, which makes the correction of WAIC unreliable, in models 1 (9), 2 (16).
#> 
#> 
#> Leave-One-Out Information Criterion (LOOIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1040   1040           977.5             1102
#>  Model 2 1018   1018           956.2             1080
#> 
#> Influential periods, whose importance sampling is unreliable, in models 1 (89), 2 (89), counted as a Pareto k above 0.


```
