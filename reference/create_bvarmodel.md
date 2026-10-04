# Create a Vector Autoregressive Model

Produces the input for the estimation of a vector autoregressive (VAR)
model.

## Usage

``` r
create_bvarmodel(
  data,
  p = 2,
  exogen = NULL,
  s = 2,
  deterministic = "const",
  seasonal = FALSE,
  iid = NULL,
  structural = FALSE,
  error = "wishart",
  quantile = 0.5,
  tvp = FALSE,
  varsel = "none",
  algorithm = NULL,
  delta_beta = 1,
  delta_sigma = 1,
  iterations = 20000,
  burnin = 2000,
  thin = 1,
  missing = "omit",
  aggregate = NULL,
  soft = NULL,
  quantile_grid = FALSE
)
```

## Arguments

- data:

  a time-series object of endogenous variables.

- p:

  an integer vector of the lag order (default is `p = 2`).

- exogen:

  an optional time-series object of external regressors.

- s:

  an integer vector of the lag order of the external regressors (default
  is `s = 2`). Ignored if `exogen` is `NULL`.

- deterministic:

  a character specifying which deterministic terms should be included.
  Available values are `"none"`, `"const"` (default) for an intercept,
  `"trend"` for a linear trend, and `"both"` for an intercept with a
  linear trend.

- seasonal:

  logical. If `TRUE`, seasonal dummy variables are generated as
  additional deterministic terms. The amount of dummies depends on the
  frequency of the time-series object provided in `data`. Defaults to
  `FALSE`.

- iid:

  an optional character vector naming endogenous variables whose
  equations carry no coefficients at all – no lags, no deterministic
  terms, nothing. They are white noise, and reach the rest of the model
  only through the error covariance. They must be the first columns of
  `data`. See 'Details'.

- structural:

  logical indicating whether data should be prepared for the estimation
  of a structural VAR model. Defaults to `FALSE`.

- error:

  character specifying the model that should be used for the estimation
  of the covariance matrix of the error term. Default is `"wishart"`.
  See 'Details'.

- quantile:

  a numeric vector of quantiles in the interval \\(0, 1)\\ that should
  be estimated. Only used, if `error = "ald"`. Defaults to `0.5`, the
  median. One model is created per quantile, so a vector produces a list
  of models, unless `quantile_grid` is `TRUE`. See 'Details'.

- tvp:

  logical indicating whether the model parameters are time varying.
  Defaults to `FALSE`.

- varsel:

  character specifying the type of variable selection algorithm that
  should be employed. Default is `"none"`. See 'Details'.

- algorithm:

  algorithm that should be used for posterior simulation. If `NULL`
  (default), the algorithm is named by `tvp` and `error`. The one
  non-standard option is `"discount"`. See 'Details'.

- delta_beta, delta_sigma:

  numeric discount factors in \\(0, 1\]\\ of the discounted model, the
  first governing the coefficients and the second the error covariance.
  Both default to 1, at which the quantity they govern does not move.
  Ignored unless `algorithm = "discount"`, and a vector in either
  produces one model per value. See 'Details'.

- iterations:

  an integer of MCMC draws excluding burn-in draws (defaults to 20000).

- burnin:

  an integer of MCMC draws used to initialize the sampler (defaults to
  2000). These draws do not enter the computation of posterior moments,
  forecasts etc.

- thin:

  an integer thinning interval of the sampler (defaults to 1). After the
  burn-in the sampler keeps the last of every `thin` draws, so it runs
  `burnin + iterations * thin` draws and still keeps `iterations`.
  Unlike
  [`thin`](https://franzmohr.github.io/bvartools/reference/thin.bvarmodel.md),
  which thins draws already made, the draws that are not kept are never
  held in memory.

- missing:

  a character, what to do with periods in which a series of `data` is
  `NA`. `"omit"` (default) drops them, as
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html) does. `"estimate"`
  keeps them and estimates what was not observed together with the
  model. See 'Details'.

- aggregate:

  an optional named list, one element per series of `data` that is
  observed at a lower frequency than the model, holding the weights with
  which it aggregates the periods of the model, oldest first, as
  [`aggregation_weights`](https://franzmohr.github.io/bvartools/reference/aggregation_weights.md)
  returns them. The series is `NA` in every period in which it is not
  observed, and its observations stand in the last period they
  aggregate. Implies `missing = "estimate"`.

- soft:

  an optional character vector of series of `data` whose observations
  hold up to a normal measurement error instead of exactly, each series
  with an error variance of its own. See 'Details'.

- quantile_grid:

  logical. If `TRUE`, the quantiles in `quantile` are estimated as one
  structural quantile VAR rather than as one model each. Needs
  `error = "ald"` and constant coefficients. Defaults to `FALSE`. See
  'Details'.

## Value

An object of class 'bvarmodel' or, if a vector is given in `p`, `s` or
`quantile`, a list of class 'modellist' with one such object per
specification. A 'bvarmodel' is a list with the elements

- `data`:

  a list with element `original`, which holds the time-series objects
  `endogen`, `exogen` and `deterministic`, and element `train`, which
  holds the estimation sample: `y`, a \\T \times K\\ time-series object
  of the endogenous variables, `x`, a time-series object of the
  regressors, and `z`, the corresponding \\TK\\ row matrix of regressors
  in SUR form, which is absent for the discounted model.

- `model`:

  a list of the specification, including `type` (`"VAR"`), `algorithm`,
  the name of the posterior simulation algorithm, `k`, `p`, `m`, `s` and
  `n`, the numbers of endogenous variables, lags, exogenous variables,
  their lags and deterministic terms, `endogen`, the names of the
  endogenous variables, and `deterministic`, `structural`, `error`,
  `varsel`, `tvp`, `iterations`, `burnin`, `thin` if it is above 1, and,
  for `error = "ald"`, `quantile` as specified.

The later steps of the workflow add the elements `priors`, `initial` and
`posterior`.

## Details

The function produces the data matrices for vector autoregressive (VAR)
models, which can also include unmodelled, non-deterministic variables:
\$\$A_0 y_t = \sum\_{i=1}^{p} A_i y\_{t - i} + \sum\_{i=0}^{s} B_i
x\_{t - i} + C d_t + u_t,\$\$ where \\y_t\\ is a K-dimensional vector of
endogenous variables, \\A_0\\ is a \\K \times K\\ coefficient matrix of
contemporaneous endogenous variables, \\A_i\\ is a \\K \times K\\
coefficient matrix of endogenous variables, \\x_t\\ is an M-dimensional
vector of exogenous regressors and \\B_i\\ its corresponding \\K \times
M\\ coefficient matrix. \\d_t\\ is an N-dimensional vector of
deterministic terms and \\C\\ its corresponding \\K \times N\\
coefficient matrix. \\p\\ is the lag order of endogenous variables,
\\s\\ is the lag order of exogenous variables, and \\u_t\\ is an error
term.

The model can be rewritten as \$\$A_0 y_t = Z_t a + u_t,\$\$ where
\\Z_t\\ is a \\KT \times K \* (Kp + M(s + 1) + N)\\ data matrix and
\\a\\ the corresponding coefficient vector. Unless structural models are
be estimated, \\A_0\\ is assumed to be an identity matrix.

If a vector is provided as argument `p` or `s`, the function will
produce a distinct model for all possible combinations of those
specifications.

If `structural = TRUE`, the data matrix \\Z_t\\ is augmented by negative
contemporaneous observations of endogenous variables, which correspond
to the coefficients in the lower triangular of \\A_0\\.

If `tvp` is `TRUE`, the respective coefficients of the above model are
assumed to be time varying. If `error` is `"sv"` or `"sv+covar"`, the
error covariance matrix is assumed to be time varying.

Argument `error` specifies the structure of the covariance matrix of the
error term and how it is estimated. Possible specifications are:

- `"wishart"`: The covariance is estimated using a Wishart prior.

- `"gamma"`: Only the diagonal elements of the covariance matrix are
  estimated using a gamma prior. Off-diagonal elements are not estimated
  and set to zero.

- `"gamma+covar"`: The diagonal elements of the covariance matrix are
  estimated using a gamma prior. Covariances are estimated based on a
  triangular decomposition.

- `"sv"`: Only the diagonal elements of the covariance matrix are
  estimated using a stochastic volatility algorithm. Off-diagonal
  elements are not estimated and set to zero.

- `"sv+covar"`: Only the diagonal elements of the covariance matrix are
  estimated using a stochastic volatility algorithm. Covariances are
  estimated based on a triangular decomposition.

- `"ald"`: The errors are assumed to follow an asymmetric Laplace
  distribution, which turns the model into a Bayesian quantile
  regression: the coefficients describe the conditional quantile
  specified in argument `quantile` instead of the conditional mean.
  Off-diagonal elements of the covariance matrix are not estimated and
  set to zero. See 'Details'.

Models with `error = "ald"` estimate a conditional quantile after Kozumi
and Kobayashi (2011). Minimising the quantile loss at \\q\\ corresponds
to maximising the likelihood of an asymmetric Laplace distribution,
which is a scale mixture of normal distributions. Conditional on the
latent scales of that mixture every equation is a weighted normal
regression, which is what makes the model a Gibbs sampler like the
others.

Three properties of these models differ from the rest of the package.
Covariances are not estimated, since rotating the equations into each
other leaves a residual whose quantile is not the one that was asked
for. A single quantile does not forecast, since the \\h\\ step ahead
quantile is not the quantile of the iterated one step ahead quantiles.
And the asymmetric Laplace is a working likelihood rather than a claim
about the data, so the posterior locates the quantile, but the spread of
the draws is not a calibrated credible interval without the adjustment
of Yang et al. (2016), which is not applied. Variable selection is
available as `"bvs"`, not as `"ssvs"`.

With `quantile_grid = TRUE` the quantiles in `quantile`, at least two,
form one model rather than one model each: the structural quantile VAR
of Chavleishvili and Manganelli (2019). Every level is estimated – each
is the chain its single-quantile model would have drawn – and the draws
of `posterior$a` and `posterior$u_scale` are stacked level by level, the
block of the first quantile first. With `structural = TRUE`, or a single
variable, the grid describes the whole distribution of each variable
given the ones ordered before it: the estimated quantiles are sorted
(Chernozhukov et al., 2010), interpolated linearly between the levels
and continued by exponential tails. That distribution is what
[`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
simulates from, variable by variable, so a grid does forecast. A
scenario given there pins variables in forecast periods, and a
`forecast_quantile` takes every draw at one level, which gives the
quantile paths the impulse responses of the model are differences of.
[`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md)
scores the density the grid describes. Only models with constant
coefficients take a grid. Summaries, plots and impulse responses work on
the single levels, which
[`split_quantile_grid`](https://franzmohr.github.io/bvartools/reference/split_quantile_grid.md)
returns.

Available specifications for argument `varsel` are:

- `"none"`: No variable selection algorithm is used.

- `"bvs"`: Bayesian variable selection as proposed in Korobilis (2013).

- `"ssvs"`: Stochastic search variable selection as proposed in George
  et al. (2008).

The one specification for argument `algorithm` is `"discount"`, the
discounted time varying parameter model of West & Harrison (1997, ch.
16) with the discounted Wishart of Uhlig (1997), estimated by
`VarTvpDiscount`. It is not a sampler: its posterior is closed form –
one pass over the sample, no chain – so `burnin` must be 0 and `thin` 1,
and `iterations` says only how many i.i.d. draws a forecast takes from
the answer. Its error covariance is the inverse Wishart whole, so
`error` must be `"wishart"` and neither variable selection nor a
structural model is available. What it buys is speed and an exact
marginal likelihood: the sum of `/posterior/loglik` is the log marginal
likelihood of the sample given the two discounts, so a grid over them
can be compared without a chain being run for any of it.

[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
estimates it like any other algorithm, the filter being part of the
vendored BayesTS core. It consumes no random numbers, so two runs agree
to the bit and a model estimated here and the same model estimated by
the `bayests` command line over a file written with
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
give the same numbers rather than merely the same distribution. What
comes back is a posterior rather than a chain: one row per period under
`posterior$a$mean`, `posterior$a$cov`, `posterior$u_sigma$scale` and
`posterior$df`, and no `coeffs` anywhere, because joining one draw per
period would look like a sampled path and is not one.

Argument `iid` restricts the equations of the variables it names to
carry no coefficients at all. Such a variable is white noise: nothing
forecasts it, it forecasts nothing, and what the model is estimated for
is its contemporaneous correlation with the errors of the equations that
do have dynamics. That is how a high-frequency surprise becomes a
variable of a monthly VAR rather than an instrument beside one,
following Jarocinski and Karadi (2020).

The restriction is exact rather than a tight prior on those
coefficients: the sampler never draws them, so they are zero in every
draw. Everything else treats the model as an ordinary VAR – the
restricted variables are still regressors in the other equations, the
error covariance still covers them, and
[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md) and
[`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md) read
the draws unchanged.

**The variables named in `iid` have to be the first columns of `data`**,
in any order among themselves. The restriction is carried by a
variable's position rather than by a list of positions, which is what
lets the sampler apply it without being told again which coefficients it
owns. A data set in another order is refused with a message naming the
columns that are there instead. It is available for models with constant
coefficients and cannot be combined with `structural` or `varsel`.

With `missing = "estimate"` the panel does not have to be observed
whole. Every observation becomes one linear constraint on the panel,
stored in `data$train$constraints`, and every sweep of the sampler draws
the periods that were not observed from their distribution given the
constraints, the coefficients and the error covariance before it draws
those (Chan, Poon and Zhu 2023). A series observed at a lower frequency
is a constraint on several periods, the weights of which `aggregate`
gives, so that a quarterly series enters a monthly model as the average,
sum or growth rate of three months (Schorfheide and Song 2015). Which
periods a constraint reaches is decided by the dates of `data`, and one
that reaches before the estimation sample is left out with a message.
The periods before the sample that the first lags reach are taken as
observed, with gaps there filled by interpolation. What `data$train$y`
holds where nothing was observed is a starting value: a linear
interpolation of what was. The draws of the completed panel come back as
`posterior$y`.

A series named in `soft` holds its constraints up to a normal error
whose precision is estimated with a gamma prior, which
[`add_prior_options`](https://franzmohr.github.io/bvartools/reference/add_prior_options.md)
sets. It suits an aggregate whose weights are an approximation, such as
the growth rate of an average.

Estimating what was not observed is available for `error = "wishart"`,
`"gamma"`, `"gamma+covar"`, `"sv"` and `"sv+covar"`, with constant or
time varying coefficients. It cannot be combined with `structural` or
the discounted model. It can be combined with `iid`, for models with
constant coefficients: a variable named there that was not observed in a
period – a high-frequency surprise whose series starts after the sample
does, as in Jarocinski and Karadi (2020) – is drawn from its white noise
given the errors of the other equations, and its gaps start at zero
rather than at an interpolation. Forecasts start from each draw's
completed panel. Forecasts are scored and conditioned on a scenario in
`data$forecast$constraints` only by models with constant coefficients
and `error = "wishart"` or `"gamma"`.

## References

Chan, J. C. C., Poon, A., & Zhu, D. (2023). High-dimensional
conditionally Gaussian state space models with missing data. *Journal of
Econometrics, 236*(1), 105468.
[doi:10.1016/j.jeconom.2023.05.005](https://doi.org/10.1016/j.jeconom.2023.05.005)

Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). *Bayesian
Econometric Methods* (2nd ed.). Cambridge: University Press.

Chavleishvili, S., & Manganelli, S. (2019). Forecasting and stress
testing with quantile vector autoregression. *ECB Working Paper*, 2330.

Chernozhukov, V., Fernandez-Val, I., & Galichon, A. (2010). Quantile and
probability curves without crossing. *Econometrica, 78*(3), 1093–1125.

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

Kozumi, H., & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
quantile regression. *Journal of Statistical Computation and Simulation,
81*(11), 1565–1578.
[doi:10.1080/00949655.2010.496117](https://doi.org/10.1080/00949655.2010.496117)

Schorfheide, F., & Song, D. (2015). Real-time forecasting with a
mixed-frequency VAR. *Journal of Business & Economic Statistics, 33*(3),
366–380.
[doi:10.1080/07350015.2014.954707](https://doi.org/10.1080/07350015.2014.954707)

Lütkepohl, H. (2006). *New Introduction to Multiple Time Series
Analysis* (2nd ed.). Berlin: Springer.

Uhlig, H. (1997). Bayesian vector autoregressions with stochastic
volatility. *Econometrica, 65*(1), 59–73.
[doi:10.2307/2171813](https://doi.org/10.2307/2171813)

West, M., & Harrison, J. (1997). *Bayesian forecasting and dynamic
models* (2nd ed.). New York: Springer.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other model set-up:
[`add_dummy_variables.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_dummy_variables.bvarmodel.md),
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_prior_options()`](https://franzmohr.github.io/bvartools/reference/add_prior_options.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`aggregation_weights()`](https://franzmohr.github.io/bvartools/reference/aggregation_weights.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 50, burnin = 10)
# Number of iterations and burnin should be much higher.
```
