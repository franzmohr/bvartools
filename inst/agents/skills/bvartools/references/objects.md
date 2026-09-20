# Model objects and posterior draws

Every `r` block on this page runs in the package's test suite, top to bottom in
one session, so the shapes asserted with `stopifnot()` are what the installed
version produces.

## The layout of a model

`?bvartools_model` is the canonical description of the object and ships with the
package. This file is the same ground in the shape an agent reads; where the two
disagree, the help page is right.

`create_bvarmodel()` and `create_bvecmodel()` return a list of class
`'bvarmodel'` or `'bvecmodel'`, and every later step returns it with something
added:

| Element | Added by | Holds |
| --- | --- | --- |
| `data$original` | `create_*model()` | The time-series objects as given: `endogen`, `exogen`, `deterministic` |
| `data$train` | `create_*model()` | The estimation sample: `y`, `x`, and `z` in SUR form; for a VEC also `w`, the lagged levels in the cointegration term |
| `model` | `create_*model()` | The specification: `k`, `p`, `m`, `s`, `n`, `endogen`, `error`, `varsel`, `tvp`, `structural`, `iterations`, `burnin`, `algorithm`; `thin` only if above 1; for a VEC also `rank`, `k_beta` |
| `priors` | `add_priors()` | See `priors.md` |
| `initial` | `add_initial_values()` | Starting values of the sampler |
| `posterior` | `add_posterior_coefficients()` and the later `add_posterior_*()` | The draws |
| `model$h`, `data$forecast` | `add_forecast_input()` | The horizon and the regressors of the forecast periods |
| `data$test$y` | `add_forecast_errors()` | What the horizon realised: one row per period and one column per variable, the periods of the horizon and nothing else. The only thing under `data` that no sampler reads |

A vector for `p`, `s`, `r` or `quantile` gives a `'modellist'` of such objects,
and `use_expanding_window()` an `'expandingwindow'`, one model per window. Every
step maps over both.

## The draws

Each block of the posterior is a list whose `coeffs` is a `coda::mcmc` matrix with
**one row per kept draw and one column per parameter**. Burn-in draws are already
gone, so there are `iterations` rows. That holds under `create_*model(thin = t)`
as well: the sampler then runs `burnin + iterations * t` draws, keeps the last of
every `t`, and labels them `t`, `2t`, ... in `coda::mcpar()`. `thin()` thins draws
that were already kept.

With `K` endogenous variables, `T` training periods and `M` coefficients,
`M = K * (K*p + m*(s + 1) + n)`, plus `K(K-1)/2` for a structural model:

| Element | Columns | Present for |
| --- | --- | --- |
| `a$coeffs` | `M`, or `T*M` with `tvp = TRUE` | every model |
| `a$lambda` | `M` | variable selection: the inclusion indicators |
| `a$sigma` | `M` | `tvp = TRUE`: the state variances |
| `u_sigma_inv$coeffs` | `K^2`, or `T*K^2` for stochastic volatility and `"ald"` | every model: the error precision matrix |
| `u_omega_inv$coeffs` | `K`, or `T*K` for stochastic volatility and `"ald"` | every error but `"wishart"`: the diagonal precisions |
| `psi$coeffs` | `K^2` | `"gamma+covar"` and `"sv+covar"`: the whole triangular matrix, not only its free elements |
| `u_scale$coeffs` | `K` | `"ald"`: the scales |
| `beta$coeffs` | `k_beta*r`, or `T*k_beta*r` with time-varying cointegration | VEC models |
| `forecast$forecasts` | `K*h` | after `add_posterior_forecasts()` |
| `forecast$errors` | `K*h` | after `add_forecast_errors()`, against `data$test$y` |
| `loglik` | `T` | after `add_posterior_loglik()` |

`loglik` is an `mcmc` matrix itself rather than a list with `coeffs`. `forecast`
is a list, but not that kind either: it is the group everything the forecast
periods produce hangs below, and its members are named after what they hold --
all of them are draws, so none of them is called `draws`. It was a matrix at
`posterior$forecast` with the errors beside it at `posterior$forecast_errors`
until bvartools 0.3.0.9000; both of those names are gone.

`posterior$loglik` is not a member of the group and is not the same statistic as
a forecast score: it evaluates each observation of the training sample under
states that have already seen it, which is what `WAIC` and `LOOIC` want.

**Never transpose a posterior draw matrix to "fix" it.** Rows are draws: a
posterior mean is `colMeans()`, a credible interval is a column quantile.
`bvar()` and `bvec()`, which collect draws from a sampler written by hand, are the
exception and take the transpose, one row per parameter.

## Reading coefficients

`a` is vec of the `K x (K*p + m*(s + 1) + n)` coefficient matrix
`[A_1 ... A_p, B_0 ... B_s, C]`, column-major:

```r
library(bvartools)
set.seed(1)

data("e1")
e1 <- diff(log(e1)) * 100
k <- ncol(e1)

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 200, burnin = 100)
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

stopifnot(identical(names(model$posterior), c("a", "u_sigma_inv")))

a <- model$posterior$a$coeffs
stopifnot(inherits(a, "mcmc"), nrow(a) == 200, ncol(a) == k * (2 * k + 1))

A <- matrix(colMeans(a), nrow = k)          # rows are equations
A1 <- A[, 1:k]                              # first lag
A2 <- A[, k + 1:k]                          # second lag
const <- A[, 2 * k + 1]
stopifnot(all(dim(A1) == c(k, k)), length(const) == k)
```

The error covariance is stored as its **inverse**. A posterior mean of the
covariance averages the inverted draws, which is not the inverse of the mean
precision:

```r
u <- model$posterior$u_sigma_inv$coeffs
stopifnot(ncol(u) == k * k)

sigma_draws <- apply(u, 1, function(x) solve(matrix(x, nrow = k)))   # k^2 by draws
Sigma <- matrix(rowMeans(sigma_draws), nrow = k)
stopifnot(all(dim(Sigma) == c(k, k)), isSymmetric(Sigma, tol = 1e-8))
```

`summary(model)` prints posterior means, standard deviations and credible bands
of every coefficient, which is usually what a report needs.

## Structural models

The contemporaneous coefficients are the last `K(K-1)/2` columns of `a`, the
free elements of the unit lower triangular `A0` stored **column by column** --
`(2,1), (3,1), ..., (K,1), (3,2), ...` -- which is the order of
`which(lower.tri(diag(K)))`. With three variables a row by row reading gives the
same matrix; from four on it does not. `A0 <- diag(K); A0[lower.tri(A0)] <- draw`
fills it correctly:

```r
structural <- create_bvarmodel(e1, p = 1, deterministic = "const",
                               structural = TRUE, error = "gamma",
                               iterations = 100, burnin = 50)
structural <- add_priors(structural, coef = list(v_i = 1),
                         sigma = list(shape = 3, rate = 1))
structural <- add_initial_values(structural)
structural <- add_posterior_coefficients(structural)

n_reduced <- k * (k + 1)                    # one lag and a constant per equation
stopifnot(ncol(structural$posterior$a$coeffs) == n_reduced + k * (k - 1) / 2,
          ncol(structural$posterior$u_omega_inv$coeffs) == k)
```

## Forecasts and the log likelihood

Forecast columns are **stacked by period**: the `K` variables of the first forecast
period, then those of the second. `predict()` reshapes them into an
`h x K x draws` array:

```r
model <- add_forecast_input(model, n_ahead = 4)
model <- add_posterior_forecasts(model)
stopifnot(model$model$h == 4)

f <- model$posterior$forecast$forecasts
stopifnot(inherits(f, "mcmc"), ncol(f) == k * 4)

pred <- predict(model, n_ahead = 4)
stopifnot(inherits(pred, "bvarprd"),
          all(dim(pred$fcst) == c(4, k, 200)),
          identical(dimnames(pred$fcst)[[2]], colnames(e1)))

by_period <- as.numeric(t(apply(pred$fcst, c(1, 2), mean)))
stopifnot(isTRUE(all.equal(unname(colMeans(f)), by_period)))

median_forecast <- apply(pred$fcst, c(1, 2), median)   # periods by variables
```

```r
model <- add_posterior_loglik(model)
tt <- nrow(model$data$train$y)
stopifnot(all(dim(model$posterior$loglik) == c(200, tt)))
```

## Thinning

`thin()` keeps every `thin`-th draw of every block:

```r
thinned <- thin(model, thin = 2)
stopifnot(nrow(thinned$posterior$a$coeffs) == 100)
```

That thins draws already made, all of which were held in memory first. To thin
in the sampler instead, give `thin` when the model is created. The chain then
runs `burnin + iterations * thin` draws and keeps the last of every `thin`, so
the posterior still has `iterations` rows and the draws that are not kept are
never stored. `coda::mcpar()` counts the draws the chain actually ran, and a
`thin()` on top counts from there:

```r
long <- create_bvarmodel(e1, p = 2, deterministic = "const",
                         iterations = 100, burnin = 100, thin = 5)
long <- add_priors(long,
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = "k", scale = 1))
long <- add_initial_values(long)
long <- add_posterior_coefficients(long)

stopifnot(long$model$thin == 5,
          nrow(long$posterior$a$coeffs) == 100,
          all(coda::mcpar(long$posterior$a$coeffs) == c(5, 500, 5)),
          all(coda::mcpar(thin(long, thin = 2)$posterior$a$coeffs) == c(10, 500, 10)))
```

Reach for it when a posterior summary moves with the seed. That means the chain
is too short for how slowly it mixes, and thinning lets it run longer without the
object growing. A chain thinned this way is exactly every `thin`-th draw of the
unthinned chain from the same seed, so if memory allows, keeping every draw is
never worse.

## VEC models

A VEC's `posterior` holds `a` (the loadings and the short-run coefficients),
`beta` and the error term. `add_posterior_forecasts()` on the `bvecmodel` adds
`forecast`, of the levels, laid out as a VAR's. Impulse responses are computed
from the VAR in levels, which `vec_to_var()` builds draw by draw; it keeps
`model`, `data` and `posterior` and drops the priors and starting values.
`recipes.md` has the whole sequence.
