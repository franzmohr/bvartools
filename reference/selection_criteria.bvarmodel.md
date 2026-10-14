# Model Selection Criteria

Calculates model selection criteria for an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
selection_criteria(object, ci = 0.95, plugin = c("mean", "best"), ...)
```

## Arguments

- object:

  an object of class 'bvarmodel'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- plugin:

  the point at which the deviance of `AIC`, `BIC` and `HQ` is evaluated.
  `"mean"`, the default, uses the posterior mean of the parameters, for
  a VEC model the best approximation of the model's rank to the
  posterior mean of \\\Pi\\. `"best"` uses the draw whose log-likelihood
  is the highest the chain visited, which is closer to the maximum these
  corrections are derived for, but grows with the length of the chain
  and is therefore comparable only across models drawn equally long.

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

## Details

`AIC`, `BIC` and `HQ` penalise the deviance at the point `plugin` names
by a count of free parameters: the coefficients of every equation, the
\\k(k + 1)/2\\ of the error covariance – or, in a structural model, its
diagonal and the free elements of \\A_0\\, which come to the same number
– for a VEC model the \\r(k + k\_{ect} - r)\\ free elements of \\\Pi =
\alpha \beta'\\, and, where the model has a volatility break, the scales
it draws and the rate their excess decays at.

**States are not counted.** The log-volatilities of a stochastic
volatility model, the coefficient paths of a time varying one and the
period scales of Student t errors are one per period, and how much of
that freedom a fit used is not a count of parameters. Those models are
compared by `WAIC` and `LOOIC`, which penalise by the flexibility the
fit actually used, or by `LPL`, which conditions on nothing a period has
not seen. A break is counted because it is a handful of parameters
rather than a path, and because the likelihood the deviance is evaluated
from scales the covariance of the periods after it.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
[`add_predictive_loglik.expandingwindow()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.expandingwindow.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`analysis_of_stored_models`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`pool_forecasts()`](https://franzmohr.github.io/bvartools/reference/pool_forecasts.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md)

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
#>  Model 1 -510.4 -509.3          -517.9             -506
#>  Model 2 -499.9 -500.5          -503.0             -496
#> 
#> 
#> Akaike Information Criterion (AIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1041   1041                                 
#>  Model 2 1031   1031                                 
#> 
#> 
#> Bayesian Information Criterion (BIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1086   1086                                 
#>  Model 2 1098   1098                                 
#> 
#> 
#> Hannan-Quinn Criterion (HQ)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1059   1059                                 
#>  Model 2 1058   1058                                 
#> 
#> 
#> Widely Applicable Information Criterion (WAIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1050   1050           982.1             1117
#>  Model 2 1042   1042           964.4             1120
#> 
#> Periods with a pointwise log-likelihood variance above 0.4, which makes the correction of WAIC unreliable, in models 1 (15), 2 (21).
#> 
#> 
#> Leave-One-Out Information Criterion (LOOIC)
#> 
#>          Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 1045   1045           979.1             1110
#>  Model 2 1028   1028           960.6             1095
#> 
#> Influential periods, whose importance sampling is unreliable, in models 1 (89), 2 (89), counted as a Pareto k above 0.


```
