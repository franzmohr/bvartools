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
`plot_forecast_errors_by_period()` plots them over the windows.

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
