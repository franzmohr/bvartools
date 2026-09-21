# Comparing models

Every `r` block on this page runs in the package's test suite, top to bottom in
one session, so what `stopifnot()` asserts is what the installed version does.

bvartools compares models in two ways: **in sample**, by criteria computed from
the pointwise log likelihood, and **out of sample**, by the errors of forecasts
from an expanding window.

## In sample

The log likelihood is not computed with the coefficients.
`add_posterior_loglik()` must come before `selection_criteria()`:

```r
library(bvartools)
set.seed(1)

data("e1")
e1 <- diff(log(e1)) * 100

models <- create_bvarmodel(e1, p = 1:3, deterministic = "const",
                           iterations = 200, burnin = 100)
stopifnot(inherits(models, "modellist"), length(models) == 3)

models <- add_priors(models,
                     coef = list(v_i = 0, v_i_det = 0),
                     sigma = list(df = "k", scale = 1))
models <- add_initial_values(models)
models <- add_posterior_coefficients(models)
models <- add_posterior_loglik(models)

# A vector of lag orders trims every model to the sample of the longest
stopifnot(length(unique(sapply(models, function(m) nrow(m$data$train$y)))) == 1)

sc <- selection_criteria(models)
stopifnot(inherits(sc, "selcritlist"), length(sc) == 3,
          all(c("LL", "AIC", "BIC", "HQ", "WAIC", "LOOIC") %in% names(sc[[1]])))
```

Each criterion is a data frame with the columns `mean`, `median`, `qlower` and
`qupper`. AIC, BIC and HQ are point estimates, so their bands are `NA`.
`choose_best_model()` returns the **position** of the best model in the list:

```r
best <- choose_best_model(sc, criterion = "WAIC")
stopifnot(best %in% 1:3)

chosen <- models[[best]]
stopifnot(inherits(chosen, "bvarmodel"))

waic <- sapply(sc, function(x) x$WAIC$mean)
stopifnot(best == which.min(waic))
```

`print(sc)` tabulates all criteria and `plot(sc, criterion = "BIC")` plots one.

Which criterion to use depends on the models:

- **AIC, BIC and HQ** count parameters, which describes a model with constant
  coefficients and a weak prior. They reproduce the textbook lag order selection.
- **WAIC and LOOIC** penalise the flexibility the fit actually used, so they are
  the ones for comparing time-varying parameters, stochastic volatility, or
  shrinkage priors such as the Minnesota prior or variable selection.
- LOOIC reports the periods whose Pareto shape is too large for its importance
  sampling to be reliable; prefer WAIC when many are flagged.

The discounted models (`algorithm = "discount"`) carry one criterion, `LML`:
the sum of `posterior$loglik`, which for them is the **exact** log marginal
likelihood of the sample rather than an estimate. There is no chain, so no
WAIC, LOOIC or information criterion; `choose_best_model()` maximises `LML`. A
vector in `delta_beta` or `delta_sigma` makes one model per value, like `p`:

```r
dm <- create_bvarmodel(e1, p = 1, deterministic = "const", algorithm = "discount",
                       delta_beta = c(0.95, 0.99, 1),
                       iterations = 50, burnin = 0, thin = 1)
dm <- add_priors(dm,
                 coef = list(v_i = 1, v_i_det = 1 / 10),
                 sigma = list(df = "k", scale = 1))
dm <- add_initial_values(dm)
dm <- add_posterior_coefficients(dm)
dm <- add_posterior_loglik(dm)

dsc <- selection_criteria(dm)
lml <- sapply(dsc, function(x) x$LML$mean)
stopifnot(length(lml) == 3, all(is.finite(lml)), is.null(dsc[[1]]$WAIC),
          choose_best_model(dsc, criterion = "LML") == which.max(lml))
```

`LML` is not comparable with `LL`, which conditions on the parameters where
`LML` integrates them out; it compares discounted models with each other.

All criteria need the models to be estimated on the **same observations**. A
vector for `p` ensures that. Models created separately, with different lag orders
or data, are put on their common sample with `combine_models()` and then
`align_model_obs()`.

## Out of sample

`use_expanding_window()` turns one specification into a series of models whose
training samples end one period apart. Call it **before** `add_priors()`, then
estimate and forecast as usual, and add the forecast errors against the full data.

`selection_criteria()` on an expanding window returns the out-of-sample criteria,
and, if the windows carry the log likelihood, the in-sample criteria of the last
window, the one estimated on the most data. Add the log likelihood as well:
some versions of the package stop without it, because the last window's forecasts
reach past the end of the data and leave it with no forecast errors.

```r
data("us_macrodata")

ew <- create_bvarmodel(us_macrodata, p = 1, deterministic = "const",
                       iterations = 100, burnin = 50)
ew <- use_expanding_window(ew, start = 2007)
stopifnot(inherits(ew, "expandingwindow"), length(ew) > 1)

ew <- add_priors(ew,
                 coef = list(v_i = 0.1, v_i_det = 0.01),
                 sigma = list(df = "k", scale = 1))
ew <- add_initial_values(ew)
ew <- add_posterior_coefficients(ew)
ew <- add_forecast_input(ew, n_ahead = 2)
ew <- add_posterior_forecasts(ew)
ew <- add_forecast_errors(ew, test_sample = us_macrodata)
ew <- add_posterior_loglik(ew)

oos <- selection_criteria(ew)
stopifnot(inherits(oos, "selcrit"),
          all(c("FE", "AFE", "RSFE", "WAIC") %in% names(oos)),
          all(c("variable", "h", "mean", "median", "qlower", "qupper") %in% names(oos$RSFE)))
```

The out-of-sample criteria are `FE`, `AFE` and `RSFE`: the forecast errors and
their absolute and root squared values, one row per variable and horizon, with
the columns `variable`, `h`, `mean`, `median`, `qlower` and `qupper`.
`plot_forecast_errors_by_period()` plots them over the windows, and
`plot(oos, criterion = "RSFE")` plots one of them by horizon, as it does for a
list of such results.

An expanding window of discounted models runs the same way. Their forecasts are
`coda::mcmc` draws like a sampler's, so `add_forecast_errors()` and
`write_to_hdf5()` take them, and `selection_criteria()` reports the forecast
error criteria beside the last window's `LML`.

`LPL`, the log predictive likelihood, is the out-of-sample criterion that is a
density rather than a distance: how likely the observations were under the
model, not how far the point forecast fell from them. It appears wherever the
densities behind it do, and they come from two places. `add_predictive_loglik()`
takes one per window of an expanding window exercise -- the criterion Koop,
Leon-Gonzalez and Strachan (2011) compare time-varying cointegration ranks with,
and the one to use where `WAIC` and `LOOIC` cannot separate models whose states
have seen the observation they are scored on. A single model carries one per
horizon of its forecast in `posterior$forecast$loglik`, written by BayesTS
against `data$test$y`; `choose_best_model(sc, criterion = "LPL")` takes the
maximum of either, `LPL` being a likelihood rather than a penalty.

Several specifications are compared by building an expanding window for each,
joining them with `combine_models()` before `add_posterior_coefficients()`, and
calling `selection_criteria()` on the combined list.
`vignette("horse-races", package = "bvartools")` runs such a comparison, and adds
forecasts published by institutions through `create_external_forecast()`.
Those have to be at the frequency of the data: annual forecasts for a quarterly
model are refused rather than read as forecasts of each year's first quarter.

## Annual forecasts against a quarterly model

Published projections -- the IMF's WEO, the ECB's, a central bank's -- are
usually **annual**: growth of annual GDP, the annual average rate of inflation,
the annual average unemployment rate. To race a quarterly model against them,
turn the model's forecasts into annual figures with `aggregate_forecasts()`
**after** `add_posterior_forecasts()`, then pass the result to
`create_external_forecast()` in place of the quarterly models.

`code` says, per variable, how the model's data were made from the untransformed
series: the seven transformation codes of FRED-MD and FRED-QD, the ones
`transform_variables()` applies. Variables left out of `code` are dropped.

| Code | Model sees | Annual figure |
| --- | --- | --- |
| 1 | `x` | a level |
| 2 | `diff(x)` | a level |
| 3 | `diff(x, differences = 2)` | a level |
| 4 | `log(x)` | a growth rate in percent |
| 5 | `diff(log(x))` | a growth rate in percent |
| 6 | `diff(log(x), differences = 2)` | a growth rate in percent |
| 7 | `diff(x / lag(x) - 1)` | a growth rate in percent |

Each draw is turned back into a path of `x`: the quarters of a year observed at
the end of a window come from the data, the rest from the draw, so every draw
becomes a draw of the annual figure. `target` picks which figure:

- `"average"` (default): growth of the annual average of `x` over the year
  before, or the annual average of a level. The WEO's convention. The growth of
  the average is also the growth of the annual sum, so GDP (a flow) and a price
  index (an average) take the same code.
- `"q4q4"`: growth of the fourth quarter over the fourth quarter of the year
  before (December over December for monthly data), or the fourth quarter's
  level. The Fed SEP's convention.

`scale` is what the output of codes 4 to 7 was multiplied by before the model
saw it: 1, the default, is `transform_variables()`'s output, 100 is log changes
in percent. Getting it wrong is silent unless `levels` is given.

`levels` is the untransformed series, `transform_variables()`'s argument `x`.
Undoing a difference needs the level it starts from, so codes **2, 3, 6 and 7
require it**; codes 1, 4 and 5 can be undone from the model's data alone, since
an unknown log level cancels from a growth rate. When given, `levels` must
reproduce the model's data under `code` and `scale` -- a wrong `scale` stops
here -- and it supplies the realised figures, as far as it reaches.

The forecast horizon must cover whole years from every origin: eight quarters
give **two** annual horizons. Horizon 1 is the year of the forecast origin, the
year of the first quarter after the training sample; horizon 2 the year after.
An external forecast is matched to the quarterly window as before -- `data_lag`
still counts quarters -- so a projection for the year of its publication is
horizon 1 whichever quarter it was published in.

```r
# The interest rate enters in first differences, inflation and unemployment as
# they are
us_data <- transform_variables(us_macrodata, c(r = 2))
us_data <- window(us_data, start = c(1959, 3))

aw <- create_bvarmodel(us_data, p = 1, deterministic = "const",
                       iterations = 100, burnin = 50)
aw <- use_expanding_window(aw, start = 2005)
aw <- add_priors(aw,
                 coef = list(v_i = 0.1, v_i_det = 0.01),
                 sigma = list(df = "k", scale = 1))
aw <- add_initial_values(aw)
aw <- add_posterior_coefficients(aw)
aw <- add_forecast_input(aw, n_ahead = 8)
aw <- add_posterior_forecasts(aw)

# Code 2 is a difference, so the levels are required
annual <- aggregate_forecasts(aw, code = c(Dp = 1, r = 2), levels = us_macrodata)
stopifnot(inherits(annual, "expandingwindow"),
          annual[[1]]$model$h == 2,
          identical(colnames(annual[[1]]$posterior$forecast$forecasts),
                    c("Dp_1", "r_1", "Dp_2", "r_2")))

# The fourth quarter instead of the annual average
q4 <- aggregate_forecasts(aw, code = c(Dp = 1, r = 2), target = "q4q4",
                          levels = us_macrodata)
stopifnot(q4[[1]]$model$aggregation$target == "q4q4")

# Annual projections for the current and the next year, one publication a
# quarter, in long format: 'period' is the year the projection is for
proj <- expand.grid(origin = 2005 + (0:7) / 4 + 0.1, ahead = 0:1,
                    variable = c("Dp", "r"), stringsAsFactors = FALSE)
proj$period <- floor(proj$origin) + proj$ahead
proj$value <- 2

ext <- create_external_forecast(proj, annual, data_lag = 1)

# No test_sample: both are scored against the annual figures of the levels
race <- add_forecast_errors(combine_models(annual, ext))
sc <- selection_criteria(race)
stopifnot(length(sc) == 2,
          identical(sc[[1]]$RSFE[, c("variable", "h")], sc[[2]]$RSFE[, c("variable", "h")]),
          all(sc[[1]]$RSFE$h %in% 1:2))
```

The realised annual values come from `levels` (without it, from the model's
data), aggregated the same way as the forecasts, and sit in `data$test$y` of
every window of both objects, so
`add_forecast_errors()` needs no `test_sample`; the windows whose years the data
do not cover yet are left without errors. To score against official annual
figures instead, pass an **annual** `ts` as `test_sample`; a quarterly one is
refused. Things that go wrong:

- Passing the quarterly models with annual forecasts stops, naming
  `aggregate_forecasts()`.
- A `modellist` with aggregated and non-aggregated models stops: aggregate all.
- Codes 2, 3, 6 or 7 without `levels` stop, and so do `levels` that do not
  reproduce the model's data -- usually `scale` 1 for data in percent.
- Aggregated models hold only forecasts and data, and a VEC comes back as a
  `bvarmodel` (its forecasts are the levels'). Aggregate last; there is nothing
  left to estimate, and `selection_criteria()` reports only `FE`, `AFE` and
  `RSFE` for them.

`vignette("macroprojections", package = "bvartools")` races models against the
published projections for Austria, aggregating the draws of `predict()` by hand.
