# The Structure of a Model Object

What the object every step of an analysis passes along is made of, where
each part comes from and how the numbers in it are laid out.

## Details

[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
and
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
return a list of class `'bvarmodel'` or `'bvecmodel'`, and every later
step returns that list with something added. Each call therefore has to
be assigned back: `model <- add_priors(model, ...)`.

The object is an ordinary list, so anything in it can be reached with
`$` or `[[`, and anything the package does not do can be done to it
directly. The layout below is the contract that makes that safe.

## The elements

- `model`:

  The specification. `k` endogenous variables, `p` lags of them, `m`
  unmodelled variables with `s` lags, `n` deterministic terms, and the
  names in `endogen`, `exogen` and `deterministic`. `error` names the
  error distribution, `varsel` the variable selection scheme, `tvp`
  whether the coefficients drift and `structural` whether the model is
  identified. `iterations` and `burnin` size the chain, `thin` is
  present only above one, and `algorithm` names the sampler that will
  draw it. A VEC model adds `rank`, `k_beta` and the deterministic terms
  restricted to the cointegration space; a quantile VAR adds `quantile`.
  [`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
  adds `seed`,
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md)
  adds `h` and `forecast_states`.

- `data$original`:

  The series as given, before any lag was taken: `endogen`, and `exogen`
  and `deterministic` where there are any. What the model is rebuilt
  from when a window of it is taken.

- `data$train`:

  The estimation sample. `y` is one row per period and one column per
  variable; `x` is the regressors in the same compact layout, one row
  per period; `z` is the same regressors in SUR form, `kron(x, I_k)`,
  which is what the samplers of a VAR read. A VEC model adds `w`, the
  lagged levels the cointegration term is formed from.

- `data$forecast`:

  Added by
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md):
  the regressors of the forecast periods, `x` in the compact layout, one
  row per horizon. A VEC model's are in the layout of its VAR
  representation, which is the form it is forecast in.

- `data$test$y`:

  Added by
  [`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md):
  what the horizon realised, one row per period and one column per
  variable, the periods of the horizon and nothing else. The only thing
  under `data` that no sampler reads – it is what a forecast is scored
  against, not anything it is estimated from – and what lets a model
  that was written to a file be scored again from the file alone.

- `priors`:

  Added by
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md):
  the prior hyperparameters, one block per parameter block of the
  posterior. Prior variances are given as precisions. Which elements are
  needed depends on `error`; the method's own help page,
  [`add_priors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md)
  or
  [`add_priors.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
  lists them.

- `initial`:

  Added by
  [`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md):
  the starting values the sampler begins from.

- `posterior`:

  Added by
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  and the later `add_posterior_*` functions. The draws; see below.

## The draws

`posterior` holds one named block per group of parameters, and each
block is a list whose `coeffs` carries the draws:

- `a`:

  The coefficients. `coeffs` has \\M\\ columns, or \\TM\\ with
  `tvp = TRUE`, the whole path of a draw stacked by period. `lambda`
  holds the inclusion indicators where variable selection is on and
  `sigma` the variances of the random walk where the coefficients drift.

- `u_sigma_inv`:

  The error precision matrix, \\K^2\\ columns per period, vectorised.
  `sigma` holds the variance of the log-volatility innovations under
  stochastic volatility, which is what a forecast steps the volatility
  forward by.

- `u_omega_inv`:

  The diagonal of that precision, \\K\\ columns per period, for every
  error but `"wishart"`.

- `psi`:

  The triangular matrix of the covariance block, \\K^2\\ columns – the
  whole matrix, not only its free elements – with `lambda` and `sigma`
  as above.

- `u_scale`:

  The scale of the asymmetric Laplace, \\K\\ columns, for a quantile
  VAR.

- `beta`:

  The cointegration vectors of a VEC model, \\k\_\beta r\\ columns, and
  `rho` where a time-varying cointegration space put a prior on its
  state autoregression.

- `q`:

  The rotations of a sign-restricted identification, added by
  [`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md).

- `loglik`:

  The pointwise log-likelihood of the estimation sample, one column per
  period, added by
  [`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md).
  A matrix of draws rather than a block with `coeffs`, since there is
  only one thing to hold.

- `forecast`:

  Everything the forecast periods produce: `forecasts`, the simulated
  paths from
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md),
  and `errors`, what
  [`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
  took against `data$test$y`. Both have \\Kh\\ columns, stacked by
  period. The members are named after what they hold rather than one of
  them being `draws`, because all of them are draws.

## Two conventions worth knowing

**Draws are rows.** Every matrix of draws is a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and one column per parameter, carrying the `start`, `end` and
`thin` of the chain. A posterior mean is
[`colMeans()`](https://rdrr.io/r/base/colSums.html) and a credible
interval a column quantile; transposing to "fix" the shape gives a
summary of the wrong thing.
[`bvar`](https://franzmohr.github.io/bvartools/reference/bvar.md) and
[`bvec`](https://franzmohr.github.io/bvartools/reference/bvec.md), which
collect the draws of a sampler written by hand, are the exception and
take the transpose.

**A period is a block of columns.** Where a quantity moves with time,
the path of one draw is stacked by period along that draw's row, so
period \\t\\ of a \\K^2\\ wide block is columns \\(t-1)K^2 + 1\\ to
\\tK^2\\. The same holds of a forecast: the \\K\\ variables of the first
horizon come first, then those of the second.

## Lists of models

Passing a vector to an argument such as `p`, `r` or `quantile` gives one
model per specification in a list of class `'modellist'`, and
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
a list of class `'expandingwindow'`, one model per window. Every step of
the workflow has a method for both and applies to each model in turn, so
the structure above is what is inside either of them.
[`open_models`](https://franzmohr.github.io/bvartools/reference/open_models.md)
is the same idea for a folder too large to hold: a handle to the files
rather than the models.

## The same layout on disk

[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
writes the object to an HDF5 file with the paths the elements above are
named by – `/data/train/y`, `/data/test/y`, `/posterior/a/coeffs`,
`/posterior/forecast/forecasts` – and
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
reads it back. It is the format the BayesTS command line runs on, which
is what makes a model estimated here and a model estimated there the
same object.

In the file every dataset is one row per quantity and one column per
draw. R's readers reverse the dimension order, so a session sees the
transpose of that, which is the draws-in-rows layout above.

## See also

[`bvartools`](https://franzmohr.github.io/bvartools/reference/bvartools-package.md)
for the workflow the steps are taken in, and
[`summary.bvarmodel`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md)
for the summaries the package takes of the draws itself.

## Examples

``` r

data("e1")
e1 <- diff(log(e1)) * 100

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
model <- add_posterior_coefficients(add_initial_values(model))

# The specification and the sample
model[["model"]][["k"]]
#> [1] 3
dim(model[["data"]][["train"]][["y"]])
#> [1] 89  3

# The draws, one row per draw and one column per parameter
draws <- model[["posterior"]][["a"]][["coeffs"]]
dim(draws)
#> [1] 20 21
colMeans(draws)
#>  [1] -0.26365467  0.04831701  0.00212573  0.35563815 -0.11413328  0.28436300
#>  [7]  0.58286268  0.25100450 -0.30652402 -0.11336451  0.05873395  0.04542951
#> [13]  0.21783670  0.02085566  0.35602898  0.34210749  0.01624599 -0.12885921
#> [19] -0.45810779  1.41104285  1.33671998
```
