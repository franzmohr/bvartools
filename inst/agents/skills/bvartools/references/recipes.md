# Worked examples

Every `r` block on this page runs in the package's test suite
(`tests/testthat/test-agent-docs.R`), top to bottom in one session. The shapes
asserted with `stopifnot()` are therefore what the installed version produces.
The draw counts are kept small so the tests stay quick. A real analysis needs
thousands of draws.

## A VAR, from data to impulse responses

```r
library(bvartools)
set.seed(42)

data("e1")
e1 <- diff(log(e1)) * 100          # log differences in percent: invest, income, cons

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 500, burnin = 200)
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),    # precisions; 0 is uninformative
                    sigma = list(df = "k", scale = 1))    # the Wishart prior on the errors
model <- add_initial_values(model)                        # from least squares by default
model <- add_posterior_coefficients(model)
```

Reading the draws. They are rows, and a parameter is a column:

```r
k <- ncol(e1)
n_x <- k * 2 + 1                   # two lags of three variables, and a constant

a <- model$posterior$a$coeffs
stopifnot(inherits(a, "mcmc"), all(dim(a) == c(500, k * n_x)))

# vec of the k x n_x matrix [A_1, A_2, c], column-major
A <- matrix(colMeans(a), nrow = k)
A1 <- A[, 1:k]                     # the first lag
stopifnot(all(dim(A1) == c(k, k)))

sigma_inv <- model$posterior$u_sigma_inv$coeffs
stopifnot(all(dim(sigma_inv) == c(500, k * k)))
```

Forecasts, then impulse responses:

```r
model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)
pred <- predict(model, n_ahead = 8)
stopifnot(inherits(pred, "bvarprd"))

oir <- irf(model, impulse = "income", response = "cons", n_ahead = 8, type = "oir")
stopifnot(inherits(oir, "bvarirf"))
```

## Comparing lag orders

A vector for `p` gives a `modellist`, and every step maps over it:

```r
models <- create_bvarmodel(e1, p = 1:3, deterministic = "const",
                           iterations = 300, burnin = 100)
stopifnot(inherits(models, "modellist"), length(models) == 3)

models <- add_priors(models,
                     coef = list(v_i = 0, v_i_det = 0),
                     sigma = list(df = "k", scale = 1))
models <- add_initial_values(models)
models <- add_posterior_coefficients(models)
models <- add_posterior_loglik(models)       # selection_criteria() needs it

crit <- selection_criteria(models)
```

`vignette("model-comparison", package = "bvartools")` covers putting models on
a common sample, and out-of-sample comparison with `use_expanding_window()`.

## A VEC, forecast in levels

```r
data("e6")
e6 <- e6 * 100

vec <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
                        iterations = 500, burnin = 200)
vec <- add_priors(vec,
                  coef = list(v_i = 0, v_i_det = 0),
                  coint = list(v_i = 0, p_tau_i = 1),    # uniform prior on the cointegration space
                  sigma = list(df = "k", scale = 1))
vec <- add_initial_values(vec)
vec <- add_posterior_coefficients(vec)

# k_beta = two variables plus the restricted constant, rank one
stopifnot(all(dim(vec$posterior$beta$coeffs) == c(500, 3)))
```

`add_forecast_input()`, `add_posterior_forecasts()` and `predict()` take the
`bvecmodel` directly. The forecast is of the levels, one row per draw and
`n_ahead * K` columns stacked by period. A VEC with time-varying coefficients or
stochastic volatility simulates its loadings, cointegration vectors and
volatility forward unless `forecast_states = "hold"`:

```r
vec <- add_forecast_input(vec, n_ahead = 8)
vec <- add_posterior_forecasts(vec)
stopifnot(all(dim(vec$posterior$forecast) == c(500, 8 * 2)))
pred <- predict(vec, n_ahead = 8)
stopifnot(inherits(pred, "bvarprd"))
```

`irf()`, `fevd()` and `spillover()` take the VAR in levels, so convert first:

```r
level_var <- vec_to_var(vec)
stopifnot(inherits(level_var, "bvarmodel"))

refused <- tryCatch(irf(vec), error = function(e) "use vec_to_var")
stopifnot(identical(refused, "use vec_to_var"))
```

## Time-varying coefficients with stochastic volatility

`tvp = TRUE` puts `shape` and `rate` for the state equation into `coef`, and
turns `v_i` and `v_i_det` into the prior precision of the coefficients before
the sample. Those must be positive: the sampler integrates that state out of the
first period's prior, which takes the inverse of its precision, so the
uninformative `v_i = 0` of a constant model is refused here. `error = "sv"`
needs the six stochastic volatility elements in `sigma`:

```r
tvp <- create_bvarmodel(e1, p = 1, deterministic = "const",
                        tvp = TRUE, error = "sv",
                        iterations = 300, burnin = 100)
tvp <- add_priors(tvp,
                  coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001),
                  sigma = list(mu = 0, v_i = 0.01, shape = 3, rate = 0.0001,
                               state_variance = 0.01, offset = 1e-5))
tvp <- add_initial_values(tvp)
tvp <- add_posterior_coefficients(tvp)

# The coefficients are a path: k(k + 1) of them in every one of tt periods,
# in period blocks
tt <- nrow(tvp$data$train$y)
stopifnot(ncol(tvp$posterior$a$coeffs) == tt * k * (k + 1))
```

## A quantile VAR

`error = "ald"` estimates a conditional quantile. It has no error covariances
and no forecast, and says so rather than producing a path:

```r
qvar <- create_bvarmodel(e1, p = 1, deterministic = "const",
                         error = "ald", quantile = 0.1,
                         iterations = 300, burnin = 100)
qvar <- add_priors(qvar,
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(shape = 3, rate = 1))
qvar <- add_initial_values(qvar)
qvar <- add_posterior_coefficients(qvar)

refused <- tryCatch(add_posterior_forecasts(add_forecast_input(qvar, n_ahead = 2)),
                    error = function(e) "refused")
stopifnot(identical(refused, "refused"))
```

## Moving a model through HDF5

`write_to_hdf5()` refuses a file that already exists. Pass `group` to put
several models into one file.

```r
path <- tempfile(fileext = ".h5")
write_to_hdf5(model, filename = path)

back <- read_model_from_hdf5(path)
stopifnot(inherits(back, "bvarmodel"),
          all(dim(back$posterior$a$coeffs) == dim(model$posterior$a$coeffs)))
```

The same file can be estimated by the BayesTS command line instead of in R.
Write the model after `add_initial_values()`, then, in a shell:

```text
bayests check model.h5
bayests posterior model.h5
```

and read the result back with `read_model_from_hdf5()`.
