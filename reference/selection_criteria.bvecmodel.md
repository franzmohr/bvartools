# Model Selection Criteria

Calculates model selection criteria for an object of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
selection_criteria(object, ci = 0.95, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- ...:

  further arguments passed to or from other methods.

## Value

A list of class 'selcrit', which also inherits the class of the model,
with the element `model` and one data frame per criterion, as described
in
[`selection_criteria.bvarmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md).

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
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data("e6")

# Create model
model <- create_bvecmodel(e6, p = 4,
                          const = "unrestricted",
                          seasonal = "unrestricted",
                          iterations = 20, burnin = 10)
#> Argument rank 'r' not specified. Generating models for r = 0, 1.
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
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
#>           Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 280.0  281.2           266.8            297.6
#>  Model 2 276.4  275.4           259.3            295.0
#> 
#> 
#> Akaike Information Criterion (AIC)
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -528.1 -528.1                                 
#>  Model 2 -515.2 -515.2                                 
#> 
#> 
#> Bayesian Information Criterion (BIC)
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -467.5 -467.5                                 
#>  Model 2 -446.7 -446.7                                 
#> 
#> 
#> Hannan-Quinn Criterion (HQ)
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -503.6 -503.6                                 
#>  Model 2 -487.5 -487.5                                 
#> 
#> 
#> Widely Applicable Information Criterion (WAIC)
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -558.8 -558.8          -560.1           -557.5
#>  Model 2 -551.1 -551.1          -552.4           -549.8
#> 
#> 
#> Leave-One-Out Information Criterion (LOOIC)
#> 
#>            Mean Median Quantile (2.5%) Quantile (97.5%)
#>  Model 1 -558.9 -558.9          -560.2           -557.7
#>  Model 2 -551.3 -551.3          -552.6           -550.0
#> 
#> Influential periods, whose importance sampling is unreliable, in models 1 (103), 2 (103), counted as a Pareto k above 0.23.


```
