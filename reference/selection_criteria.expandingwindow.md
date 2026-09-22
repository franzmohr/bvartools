# Model Selection Criteria

Calculates model selection criteria for an object of class
'expandingwindow'.

## Usage

``` r
# S3 method for class 'expandingwindow'
selection_criteria(object, ci = 0.95, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'selcrit'.

## Details

The in-sample criteria are those of the last window, the one estimated
on the most data, and the forecast error statistics those of all
windows.

If the windows hold the draws of
[`add_predictive_loglik`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
criterion `"LPL"` is the log predictive likelihood: the sum over the
windows of the log of the mean of the draws of the predictive density of
the observation that the next window adds. Its band is the normal
interval of probability `ci` around the sum with standard error
\\\sqrt{n \mathrm{Var}(lpd_t)}\\ over the \\n\\ periods, as for LOOIC.
Attribute `"terms"` holds the log predictive density of each period and
its numerical standard error, attribute `"nse"` the numerical standard
error of the sum.

## Examples

``` r

# Load data
data("e1")
orig <- diff(log(e1)) * 100
train <- window(orig, end = c(1982, 2))


# Create model
model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

model <- use_expanding_window(model, start = 1982.25)

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
#>  Variable h   MAFE RMSFE
#>      cons 1 1.1642 1.528
#>      cons 2 1.0375 1.358
#>      cons 3 0.9166 1.237
#>      cons 4    NaN   NaN
#>    income 1 1.5769 1.981
#>    income 2 1.4915 1.812
#>    income 3 1.1437 1.489
#>    income 4    NaN   NaN
#>    invest 1 3.8008 5.114
#>    invest 2 3.8017 4.518
#>    invest 3 3.5345 4.925
#>    invest 4    NaN   NaN

```
