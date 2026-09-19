# Model Selection Criteria

Calculates out-of-sample statistics for an object of class
'externalforecast'.

## Usage

``` r
# S3 method for class 'externalforecast'
selection_criteria(object, ci = 0.95, ...)
```

## Arguments

- object:

  an object of class 'externalforecast'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'selcrit'.

## Details

External forecasts are not estimated, so only out-of-sample statistics
are calculated. Since they are point forecasts, each publication
contributes a single value to the statistics of a variable and forecast
horizon.

## Examples

``` r

data("us_macrodata")

# Create model
model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "none",
                          error = "gamma", iterations = 10, burnin = 2)
# Chosen number of iterations and burn-in draws should be much higher.

model <- use_expanding_window(model, start = 2007)

# Artificial external forecasts
fcst <- expand.grid(origin = c(2007, 2007.25),
                    h = 1:2,
                    variable = c("Dp", "r"),
                    stringsAsFactors = FALSE)
fcst[["period"]] <- fcst[["origin"]] + fcst[["h"]] / 4
fcst[["value"]] <- 0

ext <- create_external_forecast(fcst, model, n_ahead = 4)

# Calculate forecast errors
ext <- add_forecast_errors(ext, test_sample = us_macrodata)

# Calculate selection criteria
selection_criteria(ext)
#> 
#> 
#> ------------------------------------------
#> Out-of-sample
#> ------------------------------------------
#>  Variable h   MAFE  RMSFE
#>        Dp 1    NaN    NaN
#>        Dp 2 0.8786 0.9128
#>        Dp 3 0.9252 0.9708
#>        Dp 4    NaN    NaN
#>         r 1    NaN    NaN
#>         r 2 5.1600 5.1608
#>         r 3 4.7850 4.7935
#>         r 4    NaN    NaN
#>         u 1    NaN    NaN
#>         u 2    NaN    NaN
#>         u 3    NaN    NaN
#>         u 4    NaN    NaN
```
