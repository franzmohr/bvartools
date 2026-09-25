---
name: bvartools
description: Write correct bvartools code — Bayesian VAR and VEC models in R with create_bvarmodel, create_bvecmodel, add_priors, add_initial_values, add_posterior_coefficients, add_forecast_input, add_posterior_forecasts, predict, irf, fevd, spillover, vec_to_var, add_posterior_loglik and selection_criteria. Use whenever R code loads bvartools, handles a 'bvarmodel', 'bvecmodel', 'modellist' or 'expandingwindow' object, or estimates a Bayesian VAR, VEC, time-varying parameter, stochastic volatility, SSVS or BVS, Minnesota prior or quantile VAR model in R. Also use when reading posterior draws, forecasts or impulse responses out of such an object, or when moving a model between R and the BayesTS command line through HDF5.
---

# Writing correct bvartools code

bvartools estimates Bayesian vector autoregressive (VAR) and vector error
correction (VEC) models. An analysis is a chain of calls on one **model object**,
each returning it with something added. The Gibbs samplers are the C++ core of
BayesTS, compiled into the package, and `set.seed()` reaches them.

```r
library(bvartools)
set.seed(42)

data("e1")
e1 <- diff(log(e1)) * 100                  # a stationary ts object, variables in columns

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 500, burnin = 200)
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)
pred <- predict(model, n_ahead = 8)
```

`iterations = 500` keeps examples quick. A real analysis needs thousands of
draws: the defaults are 20000 kept after 2000 burn-in, for a VAR and a VEC
alike. `create_*model(thin = t)` runs the chain `t` times as long and
keeps the last of every `t` draws, so a slowly mixing chain can run long while
the posterior still holds `iterations` draws; see `references/objects.md`.

## The rules that prevent wrong results

**1. Assign every step back.** `add_priors(model, ...)` on its own line computes
the priors and discards them. Always write `model <- add_priors(model, ...)`.

**2. `add_priors()` has no defaults.** `coef` and `sigma` must always be given,
`coint` for a VEC of positive rank, and `varsel` exactly when the model was created with
`varsel = "ssvs"` or `"bvs"`. Nothing is filled in: a missing required element
stops with a message naming it, and so does an element that is not recognised in
any of the four lists. `coef` needs `v_i` or a `minnesota` list. Which `sigma`
elements are required depends on `error` in `create_bvarmodel()`:

| `error` | `sigma` must contain |
| --- | --- |
| `"wishart"` (default) | `df`, `scale` |
| `"gamma"`, `"gamma+covar"` | `shape`, `rate` |
| `"sv"`, `"sv+covar"` | `mu`, `v_i`, `shape`, `rate`, `state_variance`, `offset`; with `tvp = TRUE`, `omega_v` in place of `shape` and `rate` |
| `"ald"` | `shape`, `rate` |

Prior variances are given as **precisions**: `v_i = 0` is uninformative, and a
larger `v_i` shrinks harder. `tvp = TRUE` models also need `shape` and `rate` in
`coef`, for the state equation. Afterwards, check `str(model$priors, max.level = 2)`
against what you meant.

**3. Draws are rows.** Every `model$posterior$<block>$coeffs` is a `coda::mcmc`
matrix with **one row per draw and one column per parameter**. A posterior mean
is a column mean. `a` is vec of the `K x (Kp + M(s+1) + N)` coefficient matrix
`[A_1 ... A_p, B_0 ... B_s, C]`, column-major, so `matrix(colMeans(a), nrow = K)`
rebuilds it. The same model written to HDF5 and read from Python shows the
transpose, parameters in rows, because R and HDF5 order dimensions differently.

**4. A VEC is analysed as its VAR in levels.** `irf()`, `fevd()` and
`spillover()` take a `bvarmodel`: `vec_to_var(model)` converts the estimated VEC
draw by draw, and on a `bvecmodel` those functions stop with a message saying so.
Forecasts are in levels either way. `add_forecast_input()`,
`add_posterior_forecasts()` and `predict()` also take the `bvecmodel`, and there a
model with time-varying coefficients simulates its loadings and cointegration
vectors forward, which its VAR representation cannot.

**5. Forecasts need `add_forecast_input()` first.** It sets the horizon and
builds the out-of-sample regressors. Without it, `add_posterior_forecasts()`
stops with "Model specification does not contain forecast horizon 'h'". The
`n_ahead` given there is also the longest horizon `predict()` returns. A model
with exogenous variables needs their future values in `exogen`.

**6. In-sample criteria need the log likelihood.** Call `add_posterior_loglik()`
before `selection_criteria()`, which otherwise stops.

**7. Vectors make lists of models.** A vector for `p`, `s`, `r` or `quantile`
creates one model per value in a `'modellist'`, and every step maps over it.
Index a single model with `models[[i]]`.

**8. Some combinations are refused, on purpose.** These are statistical limits,
not missing features, and the error says why:
- `structural = TRUE` with `error = "wishart"`, refused by `create_bvarmodel()`:
  `A_0` is not identified against an unrestricted covariance.
- SSVS with stochastic volatility, refused by `add_priors()`, and SSVS with
  time-varying parameters.
- `error = "ald"`, a quantile VAR, estimates no covariances and does not forecast:
  `add_posterior_forecasts()` refuses it, `irf()` takes only `type = "feir"` or
  `"custom"`, and `fevd()` and `spillover()` refuse it. Its selection scheme is
  `"bvs"` only.
- `algorithm = "discount"` needs `error = "wishart"`, `burnin = 0` and
  `thin = 1`, and takes neither variable selection nor a structural model:
  there is no chain to burn in, thin or draw an indicator along.

**9. `gen_var()` and `gen_vec()` no longer exist.** Older tutorials and answers
use them, together with `$data$Y` and `$data$SUR`. Use `create_bvarmodel()` and
`create_bvecmodel()`, whose matrices are in `model$data$train` as `y`, `x` and
`z`.

## Choosing a model

Fixed in `create_bvarmodel()` or `create_bvecmodel()`, then given priors in
`add_priors()`:

| Argument | Values |
| --- | --- |
| `error` | `"wishart"`, `"gamma"`, `"gamma+covar"`, `"sv"`, `"sv+covar"`, and `"ald"` with `quantile` for a VAR |
| `tvp` | `TRUE` for time-varying coefficients |
| `algorithm` | `NULL` (default) picks the sampler from `error` and `tvp`; `"discount"` for the closed-form discounted models (see below); `"KLGS2010"` for the VEC of Koop, Leon-Gonzalez and Strachan (2010) |
| `structural` | `TRUE` for contemporaneous coefficients (an A-model) |
| `varsel` | `"ssvs"` (George et al. 2008) or `"bvs"` (Korobilis 2013), with a `varsel` list in `add_priors()` |
| `deterministic` (VAR) | `"none"`, `"const"`, `"trend"`, `"both"` |
| `const`, `trend`, `seasonal` (VEC) | `"restricted"` to the cointegration space, or `"unrestricted"` |

`algorithm = "discount"` gives `VarTvpDiscount` or `VecTvpDiscount`: drifting
coefficients and error covariance under two discount factors, `delta_beta` and
`delta_sigma` in `(0, 1]` (one means the quantity does not move), with a
**closed-form posterior** rather than a chain. `iterations` only says how many
i.i.d. draws a forecast takes. The posterior is one row per period in
`posterior$a$mean`, `$a$cov`, `$u_sigma$scale` and `$df`, with no `coeffs`
(plain matrices, not `coda::mcmc`; the forecasts, being draws, are `mcmc`),
and the sum of `posterior$loglik` is the exact log marginal likelihood,
reported as `LML`. `coef` takes `v_i` (which must be positive), `v_i_det`,
`v_i_alpha` and `const`; a discounted VEC takes no `coint` prior, because it
conditions on a fixed cointegration matrix -- Johansen's by default, or the
`beta` argument of `add_initial_values()`.

For a VEC, `p` is the lag order **of the VAR in levels**, so the VEC itself has
`p - 1` lagged differences. `r` is the cointegration rank. The cointegration
prior goes in `coint`: `v_i` and `p_tau_i` for constant cointegration vectors,
`rho` for time-varying ones.

A Minnesota prior replaces `v_i`:
`coef = list(minnesota = list(kappa1 = 0.2, kappa2 = 0.5, kappa4 = 100))`, plus
`kappa3` when the model has exogenous variables.

## Reading the installed documentation

The help pages are the reference, and this skill only summarises them. Read the
one for the method, not the generic: `add_priors.bvarmodel` and
`add_priors.bvecmodel` document different lists. From a script or a shell:

```r
txt <- capture.output(tools::Rd2txt(
  tools::Rd_db("bvartools")[["add_priors.bvarmodel.Rd"]],
  options = list(underline_titles = FALSE)))
stopifnot(length(txt) > 0)
```

`?bvartools` has the workflow and a list of common mistakes. The vignettes
(`browseVignettes("bvartools")`) work through the Minnesota prior, SSVS,
TVP-SV VAR and VEC models, quantile VARs, sign restrictions, model comparison
and forecast horse races, including one against published annual projections.
The copy of this skill that matches the installed version is at `system.file("agents", package = "bvartools")`.

## Reference files

Read the one the task needs. Every R example in them runs in the package's test
suite, so the shapes they assert are what the installed version produces.

| File | Contents |
| --- | --- |
| `references/recipes.md` | Complete examples: a VAR, a VEC forecast in levels, a TVP-SV model, lag order comparison, a quantile VAR, a discounted VEC, and an HDF5 round trip |
| `references/priors.md` | Every element of `coef`, `sigma`, `coint` and `varsel`, what each model type requires, the Minnesota prior, and what `add_priors()` refuses |
| `references/objects.md` | The layout of a model object, the columns of every block of draws by model type, reading coefficients and covariances, forecast stacking |
| `references/analysis.md` | Forecasts with exogenous variables, the `irf()` and `fevd()` identification types, sign restrictions, spillovers |
| `references/comparison.md` | Information criteria and which to use, `choose_best_model()`, common samples, out-of-sample comparison with expanding windows |

## HDF5 and the BayesTS command line

`write_to_hdf5(model, filename, group)` writes a model in the format the
`bayests` command line reads, and `read_model_from_hdf5(filename, group)` reads
one back, with any draws the command line added. That lets a long chain run
outside R. Run `bayests check` on the file before `bayests posterior`. The file
format is documented in the BayesTS repository's `agents/` directory.

## Downstream: dfmtools

Dynamic factor models and factor augmented VARs are in the separate `dfmtools`
package. It builds on these generics and has its own skill.
