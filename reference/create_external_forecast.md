# Objects for Externally Produced Forecasts

Turns externally produced point forecasts into an object, which can be
evaluated and compared with the forecasts of Bayesian VAR and VEC
models.

## Usage

``` r
create_external_forecast(
  forecasts,
  object,
  n_ahead = NULL,
  period = "period",
  origin = "origin",
  variable = "variable",
  value = "value",
  by = NULL,
  data_lag = 1,
  select = "last"
)

# S3 method for class 'externalforecast'
add_priors(object, ...)

# S3 method for class 'externalforecast'
add_initial_values(object, ...)

# S3 method for class 'externalforecast'
add_seed(object, seed, ...)

# S3 method for class 'externalforecast'
add_posterior_coefficients(object, ...)

# S3 method for class 'externalforecast'
add_posterior_loglik(object, ...)

# S3 method for class 'externalforecast'
add_forecast_input(object, ...)

# S3 method for class 'externalforecast'
add_posterior_forecasts(object, ...)

# S3 method for class 'externalforecast'
minnesota_prior(object, ...)

# S3 method for class 'externalforecast'
ssvs_prior(object, ...)

# S3 method for class 'externalforecast'
inclusion_prior(object, ...)

# S3 method for class 'externalforecast'
thin(x, thin = 10, ...)

# S3 method for class 'externalforecast'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- forecasts:

  a data frame in long format, which contains the external point
  forecasts. See 'Details' for the required columns.

- object:

  for `create_external_forecast` a model object of class 'bvarmodel',
  'bvecmodel', 'expandingwindow' or 'modellist', which is used as the
  reference of the comparison. It provides the endogenous variables, the
  frequency of the data and the ends of the training samples, to which
  the external forecasts are matched. For the methods of the estimation
  workflow an object of class 'externalforecast'.

- n_ahead:

  the maximum forecast horizon that is considered. If `NULL` (default),
  the forecast horizon of the models in `object` is used, which requires
  that
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md)
  was already applied to them.

- period:

  name of the column of `forecasts`, which contains the period, for
  which a forecast was made.

- origin:

  name of the column of `forecasts`, which contains the period, in which
  a forecast was published.

- variable:

  name of the column of `forecasts`, which contains the names of the
  forecasted variables. They must correspond to the names of the
  endogenous variables of the models in `object`.

- value:

  name of the column of `forecasts`, which contains the values of the
  forecasts.

- by:

  name of an optional column of `forecasts`, which contains the names of
  the forecasters. If specified, one object is produced per forecaster
  and the result is a list of class 'modellist'.

- data_lag:

  an integer specifying the number of periods, by which the publication
  of the data of the endogenous variables lags behind. See 'Details'.

- select:

  either `"last"` (default) or `"first"` specifying which forecast is
  used, if multiple publications are matched to the same training
  sample. See 'Details'.

- ...:

  arguments passed forward to method.

- seed:

  not used, since external forecasts are not simulated.

- x:

  an object of class 'externalforecast'.

- thin:

  an integer specifying the thinning interval between successive draws.

- digits:

  the number of significant digits.

## Value

A list of class 'externalforecast' or, if argument `by` is specified, a
list of class 'modellist', which contains one object of class
'externalforecast' per forecaster.

## Details

Argument `forecasts` must be a data frame in long format, where each row
contains a single point forecast. The names of the required columns can
be specified in the arguments `period`, `origin`, `variable` and
`value`. Periods can be provided either as objects of class 'Date' or as
numerics, which follow the convention of
[`time`](https://rdrr.io/r/stats/time.html), i.e. 2007.25 for the second
quarter of 2007. The periods, for which a forecast was made, are rounded
to the frequency of the data of the models in `object`.

In contrast to a model, an external forecast does not have a training
sample. Therefore, each publication is matched to the training sample,
which ends closest before the publication of the forecast, so that a
model and an external forecaster are evaluated on comparable information
sets. Since the data of the endogenous variables are usually published
with a delay, argument `data_lag` can be used to specify the number of
periods, by which the last available observation lags behind the
publication of a forecast. Thus, a publication is matched to the last
training sample, which does not end after `origin - data_lag` periods.
With the default of one period an annual forecast, which was published
in the course of 2021, is matched to a training sample that ends in
2020, so that the forecast for 2021 is a one-step ahead forecast.

If multiple publications are matched to the same training sample, only
one of them is used, because otherwise a forecaster would enter the
comparison multiple times for the same period. Argument `select`
controls whether the latest (`"last"`) or the earliest (`"first"`) of
those publications is used.

The resulting object mimics the structure of a model, which was
estimated with the expanding window approach of
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md),
where each publication corresponds to one estimation window. It can be
added to a list of models with
[`combine_models`](https://franzmohr.github.io/bvartools/reference/combine_models.md)
and the functions, which are applied to obtain posterior draws, are
without effect for it. Since external forecasts are point forecasts,
their forecast errors consist of a single draw. Accordingly, the
credible bands of the out-of-sample statistics in
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
are degenerate and in-sample criteria are not available.

External forecasts do not have to be estimated, so the functions, which
add priors, initial values or posterior draws to a model are without
effect for objects of class 'externalforecast'. This allows to combine
them with model objects in a list of class 'modellist' and to apply the
usual workflow to all elements of that list.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

data("us_macrodata")

# Create model
model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "none",
                          error = "gamma", iterations = 10, burnin = 2)
# Chosen number of iterations and burn-in draws should be much higher.

model <- use_expanding_window(model, start = 2007)

# Artificial external forecasts of two forecasters
fcst <- expand.grid(origin = c(2007, 2007.25),
                    h = 1:2,
                    variable = c("Dp", "r"),
                    forecaster = c("A", "B"),
                    stringsAsFactors = FALSE)
fcst[["period"]] <- fcst[["origin"]] + fcst[["h"]] / 4
fcst[["value"]] <- 0

# Create objects of the external forecasts, where the data of the endogenous
# variables are assumed to be published with a delay of one quarter
ext <- create_external_forecast(fcst, model, n_ahead = 4,
                                by = "forecaster", data_lag = 1)

# Calculate forecast errors
ext <- add_forecast_errors(ext, test_sample = us_macrodata)

# Compare the forecast performance
selection_criteria(ext)
#> 
#> 
#> ------------------------------------------
#> Out-of-sample
#> ------------------------------------------
#> 
#> Mean absolute forecast errors (MAFE)
#> 
#>  Variable h Model 1 Model 2
#>        Dp 1     NaN     NaN
#>         u 1     NaN     NaN
#>         r 1     NaN     NaN
#>        Dp 2  0.8786  0.8786
#>         u 2     NaN     NaN
#>         r 2  5.1600  5.1600
#>        Dp 3  0.9252  0.9252
#>         u 3     NaN     NaN
#>         r 3  4.7850  4.7850
#>        Dp 4     NaN     NaN
#>         u 4     NaN     NaN
#>         r 4     NaN     NaN
#> 
#> 
#> Root mean squared forecast errors (RMSFE)
#> 
#>  Variable h Model 1 Model 2
#>        Dp 1     NaN     NaN
#>         u 1     NaN     NaN
#>         r 1     NaN     NaN
#>        Dp 2  0.9128  0.9128
#>         u 2     NaN     NaN
#>         r 2  5.1608  5.1608
#>        Dp 3  0.9708  0.9708
#>         u 3     NaN     NaN
#>         r 3  4.7935  4.7935
#>        Dp 4     NaN     NaN
#>         u 4     NaN     NaN
#>         r 4     NaN     NaN
```
