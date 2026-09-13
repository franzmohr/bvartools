# Forecasts, impulse responses, variance decompositions and spillovers

Every `r` block on this page runs in the package's test suite, top to bottom in
one session, so the shapes asserted with `stopifnot()` are what the installed
version produces.

All four take a `'bvarmodel'` with posterior draws. **Convert a VEC with
`vec_to_var()` first.** For a model with time-varying parameters or stochastic
volatility, `period` picks the period whose draws are used; the default is the
last one.

```r
library(bvartools)
set.seed(1)

data("us_macrodata")                         # Dp (inflation), u (unemployment), r (Fed funds rate)
k <- ncol(us_macrodata)

model <- create_bvarmodel(us_macrodata, p = 2, deterministic = "const",
                          iterations = 300, burnin = 100)
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 0.1),
                    sigma = list(df = "k", scale = 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
```

## Forecasts

`add_forecast_input()` fixes the horizon and builds the regressors of the
forecast periods; `add_posterior_forecasts()` simulates them; `predict()`
summarises. Deterministic terms are extended automatically. Future values of
**exogenous variables** must be passed to `add_forecast_input()` in `exogen`.

```r
model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)
pred <- predict(model, n_ahead = 8)

stopifnot(inherits(pred, "bvarprd"), all(dim(pred$fcst) == c(8, k, 300)))
```

`pred$fcst` holds the draws, periods by variables by draws, and `pred$y` the
training data; `plot(pred, n_pre = 20)` shows the last 20 observations with the
bands. `n_ahead` in `predict()` cannot exceed the one given to
`add_forecast_input()`; a longer request is cut back to it.

## Impulse responses

`irf(model, impulse, response, n_ahead, type)` returns a time series from period 0
to `n_ahead`, with the lower bound, median and upper bound of the credible band
named after their quantiles:

| `type` | Identification |
| --- | --- |
| `"feir"` (default) | Forecast error: a unit shock to the reduced-form error |
| `"oir"` | Orthogonalised, by the Cholesky factor: the ordering of the variables matters |
| `"gir"` | Generalised (Pesaran and Shin 1998): order-invariant |
| `"sign"` | Sign restricted, after `add_sign_restrictions()` |
| `"custom"` | An impact matrix given in `impact` |
| `"sir"`, `"sgir"` | Structural and structural generalised, the only two for `structural = TRUE` |

```r
oir <- irf(model, impulse = "r", response = "Dp", n_ahead = 8, type = "oir")
stopifnot(inherits(oir, "bvarirf"), nrow(oir) == 9,
          identical(colnames(oir), c("2.5%", "50%", "97.5%")))

draws <- irf(model, impulse = "r", response = "Dp", n_ahead = 8,
             type = "oir", keep_draws = TRUE)
stopifnot(all(dim(draws) == c(300, 9)))     # draws by periods 0 to 8
```

`ci` sets the band (default `0.95`), `shock` its size and `cumulative = TRUE`
accumulates the responses.

## Variance decompositions

`fevd(model, response, n_ahead, type)` returns one column per shock and one row
per period from 0. Its default type is `"oir"`, whose shares sum to one; `"gir"`
shares only do with `normalise_gir = TRUE`.

```r
vd <- fevd(model, response = "Dp", n_ahead = 8)
stopifnot(inherits(vd, "bvarfevd"),
          identical(colnames(vd), colnames(us_macrodata)),
          all(abs(rowSums(vd) - 1) < 1e-6))
```

## Sign restrictions

A data frame with one row per restriction: the `impulse` variable naming the
shock, the `response`, the `sign` (`1` or `-1`) and the `horizon`.
`add_sign_restrictions()` searches rotations for each draw and stores the accepted
ones; `irf()` and `fevd()` then use them with `type = "sign"`.

```r
restrictions <- data.frame(
  impulse  = "r",
  response = c("r", "Dp", "u"),
  sign     = c(1, -1, 1),
  horizon  = rep(0:1, each = 3)
)

model <- add_sign_restrictions(model, restrictions, max_tries = 500)
sign_ir <- irf(model, impulse = "r", response = "Dp", n_ahead = 8, type = "sign")
stopifnot(inherits(sign_ir, "bvarirf"), nrow(sign_ir) == 9)
```

Draws for which no admissible rotation is found within `max_tries` are dropped
from the result; a restriction set that admits almost nothing leaves few draws,
so check `model$model$sign_restrictions` before relying on the bands.

## Spillovers

`spillover()` computes the connectedness measures of Diebold and Yilmaz (2012)
from generalised variance decompositions by default:

```r
sp <- spillover(model, n_ahead = 8)
stopifnot(inherits(sp, "bvarspillover"),
          all(c("total", "to", "from", "net", "table") %in% names(sp)))
```

A structural model is not supported, since the measures decompose the reduced
form.
