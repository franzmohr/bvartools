# Adaptive Priors, Stationarity and the Steady State

Adds the options of the prior that act on top of the one
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
sets: hyperparameters estimated from the data, a restriction of the
coefficients to a stationary model, a prior on the unconditional mean in
place of one on the intercept, and the prior of the error of soft
constraints.

## Usage

``` r
add_prior_options(
  object,
  shrinkage = "none",
  stationary = FALSE,
  steady_state = NULL,
  constraints = NULL
)
```

## Arguments

- object:

  an object of class 'bvarmodel' or 'modellist', usually the result of
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

- shrinkage:

  a character, `"none"` (default), `"minnesota"`, `"horseshoe"` or
  `"normal_gamma"`, or a list with element `type` naming one of them and
  elements `group`, `shape`, `rate` and, for `"normal_gamma"`, `theta`
  or `theta_rate`. See 'Details'.

- stationary:

  logical. If `TRUE`, every draw of the coefficients is a stationary
  model. Defaults to `FALSE`.

- steady_state:

  an optional list with elements `mu`, the prior mean of the
  unconditional mean of the endogenous variables, and `v_i`, its prior
  precision, each a number or one value per variable. See 'Details'.

- constraints:

  an optional list with elements `shape` and `rate`, the gamma prior of
  the error precision of soft constraints, a number or one value per
  series named in argument `soft` of
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).
  Defaults to `list(shape = 3, rate = 0.01)` where the model has soft
  constraints.

## Value

The object in `object` with `model$shrinkage`, `model$stationary` and
`model$steady_state` set where they apply, and the priors in
`priors$a$shrinkage`, `priors$mu` and `priors$constraints`.

## Details

Argument `shrinkage` lets the data decide how tightly the prior of
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
shrinks the coefficients of the lagged endogenous variables. Their prior
variances are multiplied by scales that are drawn with the coefficients:

- `"minnesota"`: one scale per group of coefficients with an inverse
  gamma prior of `shape` and `rate`, the hierarchical Minnesota prior of
  Chan (2021). By default own lags form group 1 and the lags of the
  other variables group 2, with `shape = 3` and `rate = 2`, a prior mean
  of one: the prior of `add_priors` is the centre of the one estimated.

- `"horseshoe"`: the horseshoe prior of Carvalho, Polson and Scott
  (2010), a half-Cauchy scale per group and one per coefficient, drawn
  as in Makalic and Schmidt (2016). By default all lags form one group.

- `"normal_gamma"`: the normal-gamma prior of Griffin and Brown (2010)
  as Huber and Feldkircher (2019) put it on a VAR. Each coefficient's
  scale \\\psi_j\\ has a gamma prior of shape `theta` and rate \\\theta
  \lambda_g / 2\\, and each group's \\\lambda_g\\ a gamma prior of
  `shape` and `rate`. By default the lags of order \\l\\ form group
  \\l\\, with `theta = 0.1` and `shape = rate = 0.01`. `theta` is held
  fixed; the smaller it is, the more the prior pushes small coefficients
  to zero while leaving large ones alone. Give `theta_rate` instead to
  draw \\\theta\\ under an exponential prior of that rate, by a
  Metropolis-Hastings step tuned during the burn-in, as Huber and
  Feldkircher (2019) do; the draws come back as `posterior$a$theta`, and
  the chain needs a burn-in to tune the step.

Coefficients in group 0 – by default the deterministic terms and the
exogenous variables – keep the prior of `add_priors`. A `group` of one's
own has one element per coefficient, in the order of the columns of
`data$train$z`. The prior must be diagonal wherever it shrinks.

Argument `stationary` redraws a draw of the coefficients that is not
stationary up to 100 times, and keeps the previous draw if none of them
is, which leaves the posterior restricted to the stationary region
invariant.

Argument `steady_state` puts the prior on the unconditional mean \\\mu\\
of the endogenous variables rather than on the intercept, which is
\\(I - \sum_i A_i) \mu\\ in every draw, following Villani (2009). A
prior belief about the level a series returns to is more often available
than one about an intercept. The model must have lags and an intercept
and nothing else, and cannot be combined with `shrinkage`. The draws of
\\\mu\\ come back as `posterior$mu`.

The three are available for models with constant coefficients and
`error = "wishart"`, `"gamma"`, `"gamma+covar"`, `"sv"` or `"sv+covar"`,
without variable selection or a structural form.

## References

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
estimator for sparse signals. *Biometrika, 97*(2), 465–480.
[doi:10.1093/biomet/asq017](https://doi.org/10.1093/biomet/asq017)

Griffin, J. E., & Brown, P. J. (2010). Inference with normal-gamma prior
distributions in regression problems. *Bayesian Analysis, 5*(1),
171–188. [doi:10.1214/10-BA507](https://doi.org/10.1214/10-BA507)

Huber, F., & Feldkircher, M. (2019). Adaptive shrinkage in Bayesian
vector autoregressive models. *Journal of Business & Economic
Statistics, 37*(1), 27–39.
[doi:10.1080/07350015.2016.1256217](https://doi.org/10.1080/07350015.2016.1256217)

Chan, J. C. C. (2021). Minnesota-type adaptive hierarchical priors for
large Bayesian VARs. *International Journal of Forecasting, 37*(3),
1212–1226.
[doi:10.1016/j.ijforecast.2021.01.002](https://doi.org/10.1016/j.ijforecast.2021.01.002)

Makalic, E., & Schmidt, D. F. (2016). A simple sampler for the horseshoe
estimator. *IEEE Signal Processing Letters, 23*(1), 179–182.
[doi:10.1109/LSP.2015.2503725](https://doi.org/10.1109/LSP.2015.2503725)

Villani, M. (2009). Steady-state priors for vector autoregressions.
*Journal of Applied Econometrics, 24*(4), 630–650.
[doi:10.1002/jae.1065](https://doi.org/10.1002/jae.1065)

## See also

Other model set-up:
[`add_dummy_variables.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_dummy_variables.bvarmodel.md),
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`aggregation_weights()`](https://franzmohr.github.io/bvartools/reference/aggregation_weights.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r
data("e1")
e1 <- diff(log(e1)) * 100

model <- create_bvarmodel(e1, p = 2, iterations = 20, burnin = 10)
model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))
model <- add_prior_options(model, shrinkage = "minnesota", stationary = TRUE)
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
```
