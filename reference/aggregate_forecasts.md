# Aggregate Forecasts to Annual Figures

Turns the forecast draws of models estimated on quarterly or monthly
data into draws of the annual figures they imply, so that the models can
be compared with annual forecasts, such as the projections of
international institutions, on the same annual basis.

## Usage

``` r
aggregate_forecasts(object, code, target = "average", levels = NULL, scale = 1)
```

## Arguments

- object:

  an object of class 'bvarmodel', 'bvecmodel', 'expandingwindow' or
  'modellist', whose models contain forecasts, i.e.
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
  was already applied to them.

- code:

  a named integer vector of the transformation codes of FRED-MD and
  FRED-QD, `1:7`, which were applied to the variables of the models –
  see
  [`transform_variables`](https://franzmohr.github.io/bvartools/reference/transform_variables.md)
  and 'Details'. The names are endogenous variables of the models.
  Endogenous variables, which are not named, are dropped from the
  comparison.

- target:

  either `"average"` (default) or `"q4q4"`, the annual figure that is
  compared. See 'Details'.

- levels:

  an optional time-series object with the untransformed series, to which
  `code` was applied, i.e. the argument `x` of
  [`transform_variables`](https://franzmohr.github.io/bvartools/reference/transform_variables.md).
  Its columns are named after the variables. Required for codes 2, 3, 6
  and 7. See 'Details'.

- scale:

  the factor, by which the transformed series of codes 4 to 7 were
  multiplied before they were passed to the models, such as 100 for log
  differences in percent. The default of 1 corresponds to the result of
  [`transform_variables`](https://franzmohr.github.io/bvartools/reference/transform_variables.md).

## Value

An object of the same class as `object`, whose models contain the draws
of the annual forecasts in `posterior$forecast$forecasts`.

## Details

Every draw of a forecast is turned back into a path of the untransformed
series by reversing the transformation of its code, and continued from
the periods, which were observed at the end of the training sample. The
periods of a year, which were already observed, are thus taken from the
data and the remaining periods from the draw, so that every draw of the
forecast becomes a draw of the annual figure.

The codes determine the annual figure. The multiplicative codes 4 to 7 –
logarithms, their differences and differences of growth rates – describe
a series such as GDP or a price index, whose annual figure is a growth
rate in percent. The additive codes 1 to 3 describe a series such as an
unemployment rate or an interest rate, whose annual figure is a level.
Argument `target` determines, which growth rate and which level:

- `"average"`:

  the growth of the annual average of the levels over the annual average
  of the year before for codes 4 to 7, and the annual average for codes
  1 to 3. This is the convention of, e.g., the World Economic Outlook of
  the IMF. The growth of the annual average is also the growth of the
  annual sum, since both years have the same number of periods, so it
  serves flows such as GDP and averages such as prices alike.

- `"q4q4"`:

  the growth of the level of the last period of the year over the last
  period of the year before, i.e. the fourth quarter over the fourth
  quarter or December over December, for codes 4 to 7, and the level of
  the last period of the year for codes 1 to 3. This is the convention
  of, e.g., the Summary of Economic Projections of the Federal Reserve.

Reversing a difference requires the level, from which it starts. For
codes 2, 3, 6 and 7 it is taken from argument `levels`, which is
therefore required for them. Codes 1, 4 and 5 can be reversed from the
data of the models, because the unknown level of a logarithm cancels
from a growth rate. If `levels` is provided, it must reproduce the data
of the models under `code` and `scale`, which guards against a wrong
`scale`, and it also provides the realised annual figures.

Horizon 1 is the year of the forecast origin, i.e. the year of the first
period after the training sample, and horizon 2 the year after it. The
number of annual horizons is the number of years that the forecast
horizon covers from every origin, `floor(n_ahead / frequency)`, so that
a quarterly model, which forecasts eight quarters, provides the current
and the next year.

The realised annual figures, against which the aggregated forecasts are
scored, are obtained from `levels` or, without it, from the data of the
models in the same way, and are put in `data$test$y`.
[`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
uses them, if its argument `test_sample` is omitted. Otherwise,
`test_sample` must be an annual time series, for example of official
annual figures.

The result can be passed as argument `object` to
[`create_external_forecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
which then reads the periods of the external forecasts as years, but
matches their publications to the training samples of the models at the
frequency of the data. The external forecasts are scored against the
same realised annual figures as the models.

The aggregated models only contain the forecasts and the data, which are
required to evaluate them. Therefore, aggregation is the last step
before
[`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
and
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md),
which only provides out-of-sample statistics for them. Aggregated VEC
models become objects of class 'bvarmodel', because their forecasts are
those of the levels.

## References

McCracken, M. W., & Ng, S. (2021). FRED-QD: A quarterly database for
macroeconomic research. *Federal Reserve Bank of St. Louis Review,
103*(1), 1–44.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
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
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

data("us_macrodata")

# Inflation and the interest rate enter the model as they are
model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "const",
                          iterations = 10, burnin = 2)
# Chosen number of iterations and burn-in draws should be much higher.

model <- use_expanding_window(model, start = 2005)
model <- add_priors(model, coef = list(v_i = 0.1, v_i_det = 0.01),
                    sigma = list(df = "k", scale = 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)

# Annual averages of inflation and the interest rate
annual <- aggregate_forecasts(model, code = c(Dp = 1, r = 1))

# Artificial annual forecasts of the current and the next year, published in
# the middle of each quarter
fcst <- expand.grid(origin = 2005 + (0:7) / 4 + 0.1, year = 0:1,
                    variable = c("Dp", "r"), stringsAsFactors = FALSE)
fcst[["period"]] <- floor(fcst[["origin"]]) + fcst[["year"]]
fcst[["value"]] <- 2

ext <- create_external_forecast(fcst, annual, data_lag = 1)

# Both are scored against the annual averages of the data
models <- add_forecast_errors(combine_models(annual, ext))
selection_criteria(models)
#> 
#> 
#> ------------------------------------------
#> Out-of-sample
#> ------------------------------------------
#> 
#> Mean absolute forecast errors (MAFE)
#> 
#>  Variable h Model 1 Model 2
#>        Dp 1  0.3278   1.306
#>         r 1  0.7703   2.090
#>        Dp 2  0.4959   1.263
#>         r 2  1.7598   2.994
#> 
#> 
#> Root mean squared forecast errors (RMSFE)
#> 
#>  Variable h Model 1 Model 2
#>        Dp 1  0.4274   1.322
#>         r 1  1.0754   2.267
#>        Dp 2  0.6196   1.287
#>         r 2  2.1979   2.994
```
