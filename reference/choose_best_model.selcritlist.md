# Choose Best Model

Chooses the best model according the selection criteria in an object of
class 'selcritlist'.

## Usage

``` r
# S3 method for class 'selcritlist'
choose_best_model(object, criterion = "WAIC", ...)
```

## Arguments

- object:

  object of class 'selcritlist', usually, a result of a call to
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md).

- criterion:

  the selection criterion that should be used. Available choices are
  `"LL"`, `"AIC"`, `"BIC"`, `"HQ"`, `"WAIC"` (default), `"LOOIC"` and
  `"LPL"`.

- ...:

  further arguments passed to or from other methods.

## Value

An integer giving the position of the best model in the list provided in
argument `object`.

## Details

If argument `criterion` is "LL" or "LPL", the model with the maximum
value is chosen, otherwise, the model with the minimum value.

Which criterion to use depends on the models that are compared. See
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md),
whose details set out what each of them penalises and when a count of
parameters ceases to describe a model. The default is `"WAIC"`, because
it is the criterion that stays defined across the models this package
estimates: it penalises by the flexibility the fit used rather than by a
count of parameters, which describes neither a model whose coefficients
or variances follow a state equation nor one whose prior shrinks them.
`"AIC"`, `"BIC"` and `"HQ"` remain the criteria for the lag order of a
model with constant coefficients and a weak prior, which is the case
they are derived for. `"LPL"`, the log predictive likelihood of an
expanding window exercise with
[`add_predictive_loglik`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
is the criterion for the rank of a VEC model whose coefficients or
variances follow a state equation.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r
# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 1:3, deterministic = "const",
                          iterations = 50, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Obtain posterior log-likelihoods
model <- add_posterior_loglik(model)

# Calculate selection criteria
sel <- selection_criteria(model)

# View results
sel
#> 
#> ------------------------------------------
#> In-sample
#> ------------------------------------------
#> 
#> Log-likelihood
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -506.1 -506.2          -511.8           -501.2
#>  Model 2 -496.7 -496.3          -506.4           -491.2
#>  Model 3 -497.2 -497.1          -506.0           -489.7
#> 
#> 
#> Akaike Information Criterion (AIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1029   1029                                 
#>  Model 2 1021   1021                                 
#>  Model 3 1030   1030                                 
#> 
#> 
#> Bayesian Information Criterion (BIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1074   1074                                 
#>  Model 2 1088   1088                                 
#>  Model 3 1120   1120                                 
#> 
#> 
#> Hannan-Quinn Criterion (HQ)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1047   1047                                 
#>  Model 2 1048   1048                                 
#>  Model 3 1066   1066                                 
#> 
#> 
#> Widely Applicable Information Criterion (WAIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1037   1037           968.4             1106
#>  Model 2 1025   1025           960.8             1089
#>  Model 3 1038   1038           973.3             1102
#> 
#> Periods with a pointwise log-likelihood variance above 0.4, which makes the correction of WAIC unreliable, in models 1 (14), 2 (21), 3 (31).
#> 
#> 
#> Leave-One-Out Information Criterion (LOOIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1035   1035           967.4             1102
#>  Model 2 1022   1022           959.1             1085
#>  Model 3 1033   1033           971.1             1095
#> 
#> Influential periods, whose importance sampling is unreliable, in models 1 (16), 2 (31), 3 (32), counted as a Pareto k above 0.41.

# Choose best model according to WAIC, the default
choose_best_model(sel)
#> [1] 2

# Choose best model according to AIC
choose_best_model(sel, criterion = "AIC")
#> [1] 2
```
