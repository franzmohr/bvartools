# Create Vector Error Correction Models

Produces the input for the estimation of a vector error correction (VEC)
model.

## Usage

``` r
create_bvecmodel(
  data,
  p,
  exogen = NULL,
  s = NULL,
  r = NULL,
  const = NULL,
  trend = NULL,
  seasonal = NULL,
  structural = FALSE,
  error = "wishart",
  tvp = FALSE,
  varsel = "none",
  algorithm = NULL,
  delta_beta = 1,
  delta_sigma = 1,
  iterations = 20000,
  burnin = 2000,
  thin = 1
)
```

## Arguments

- data:

  a time-series object of endogenous variables.

- p:

  an integer vector of the lag order of the endogenous variables in the
  corresponding VAR in levels, which must be at least 1. Following
  common convention, the resulting VEC model contains \\p - 1\\ lags of
  the differenced endogenous variables. There is no default. See
  'Details'.

- exogen:

  an optional time-series object of external regressors.

- s:

  an integer vector of the lag order of the exogenous variables in the
  corresponding VAR in levels, which must be at least 1. The resulting
  VEC model contains the contemporaneous difference of the exogenous
  variables and \\s - 1\\ lags of it. Must be specified if `exogen` is
  given and is ignored otherwise. See 'Details'.

- r:

  an integer vector of the cointegration rank. See 'Details'.

- const:

  a character specifying whether a constant term enters the error
  correction term (`"restricted"`) or the non-cointegration term as an
  `"unrestricted"` variable. If `NULL` (default) no constant term will
  be added.

- trend:

  a character specifying whether a trend term enters the error
  correction term (`"restricted"`) or the non-cointegration term as an
  `"unrestricted"` variable. If `NULL` (default) no constant term will
  be added.

- seasonal:

  a character specifying whether seasonal dummies should be included in
  the error correction term (`"restricted"`) or in the
  non-cointegreation term as `"unrestricted"` variables. If `NULL`
  (default) no seasonal terms will be added. The amount of dummy
  variables will be automatically detected and depends on the frequency
  of the time-series object provided in `data`.

- structural:

  logical indicating whether data should be prepared for the estimation
  of a structural VEC model, whose lower triangular matrix \\A_0\\
  multiplies the differences of the endogenous variables.

- error:

  character specifying the model that should be used for the estimation
  of the covariance matrix of the error term. Default is `"wishart"`.
  See 'Details'.

- tvp:

  logical indicating whether the model parameters are time varying.

- varsel:

  character specifying the type of variable selection algorithm that
  should be employed. Default is `"none"`. See 'Details'.

- algorithm:

  algorithm that should be used for posterior simulation. If `NULL`
  (default), standard algorithms will be used. See 'Details' for
  available non-standard options.

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
  [`thin`](https://franzmohr.github.io/bvartools/reference/thin.bvecmodel.md),
  which thins draws already made, the draws that are not kept are never
  held in memory.

## Value

An object of class 'bvecmodel' or, if a vector is given in `p`, `s` or
`r`, a list of class 'modellist' with one such object per specification.
A 'bvecmodel' is a list with the elements

- `data`:

  a list with element `original`, which holds the time-series objects
  `endogen`, `exogen` and `deterministic` in levels, and element
  `train`, which holds the estimation sample: `y`, a \\T \times K\\
  time-series object of the differenced endogenous variables, `w`, the
  lagged levels that enter the cointegration term, `x`, the remaining
  regressors, and `z`, the corresponding \\TK\\ row matrix of regressors
  in SUR form, which is absent for the discounted model.

- `model`:

  a list of the specification with the same elements as that of a
  'bvarmodel', where `type` is `"VEC"` and `p` is the lag order of the
  VAR in levels, and additionally `rank`, the cointegration rank,
  `k_beta`, the number of variables in the cointegration term, and
  `n_restricted`, the number of deterministic terms restricted to it.

The later steps of the workflow add the elements `priors`, `initial` and
`posterior`.

## Details

The function produces the variable matrices of vector error correction
(VEC) models, which can also include exogenous variables: \$\$\Delta y_t
= \Pi w_t + \sum\_{i=1}^{p-1} \Gamma\_{i} \Delta y\_{t - i} +
\sum\_{i=0}^{s-1} \Upsilon\_{i} \Delta x\_{t - i} + C^{UR} d^{UR}\_t +
u_t,\$\$ where \\\Delta y_t\\ is a \\K \times 1\\ vector of differenced
endogenous variables, \\w_t\\ is a \\(K + M + N^{R}) \times 1\\ vector
of cointegration variables, \\\Pi\\ is a \\K \times (K + M + N^{R})\\
matrix of cointegration parameters, \\\Gamma_i\\ is a \\K \times K\\
coefficient matrix of endogenous variables, \\\Delta x_t\\ is a \\M
\times 1\\ vector of differenced exogenous regressors, \\\Upsilon_i\\ is
a \\K \times M\\ coefficient matrix of exogenous regressors,
\\d^{UR}\_t\\ is a \\N^{UR} \times 1\\ vector of deterministic terms,
and \\C^{UR}\\ is a \\K \times N^{UR}\\ coefficient matrix of
deterministic terms that do not enter the cointegration term. \\p\\ is
the lag order of endogenous variables and \\s\\ is the lag order of
exogenous variables of the corresponding VAR model. \\u_t\\ is a \\K
\times 1\\ error term.

If an integer vector is provided as argument `p`, `s` or `r`, the
function will produce a distinct model for all possible combinations of
those specifications.

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

Available specifications for argument `varsel` are:

- `"none"`: No variable selection algorithm is used.

- `"bvs"`: Bayesian variable selection as proposed in Korobilis (2013).

- `"ssvs"`: Stochastic search variable selection as proposed in George
  et al. (2008).

Available specifications for argument `algorithm` are:

- `"KLGS2010"`: Algorithm proposed in Koop, León-González & Strachan
  (2010).

- `"discount"`: The discounted time varying parameter model of West &
  Harrison (1997, ch. 16) with the discounted Wishart of Uhlig (1997),
  estimated by `VecTvpDiscount`.

The discounted model is not a sampler. Its posterior is closed form –
one pass over the sample, no chain – so `burnin` must be 0 and `thin` 1,
and `iterations` says only how many i.i.d. draws a forecast takes from
the answer. Its error covariance is the inverse Wishart whole, so
`error` must be `"wishart"` and neither variable selection nor a
structural model is available. What it buys is speed and an exact
marginal likelihood: the sum of `/posterior/loglik` is the log marginal
likelihood of the sample given the rank, the cointegration matrix and
the two discounts, so a grid over any of them can be compared without a
chain being run for any of it.

The one assumption that separates it from the sampling VEC models is
that the cointegration space is given rather than estimated.
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
puts Johansen's maximum likelihood estimate at `/initial/beta` and the
model conditions on it, so what drifts is the adjustment to the long-run
relations and not the relations themselves. That is a different question
from the one `"KLGS2010"` answers, not a cheaper way of answering the
same one.

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

## References

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224–242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Uhlig, H. (1997). Bayesian vector autoregressions with stochastic
volatility. *Econometrica, 65*(1), 59–73.
[doi:10.2307/2171813](https://doi.org/10.2307/2171813)

West, M., & Harrison, J. (1997). *Bayesian forecasting and dynamic
models* (2nd ed.). New York: Springer.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

# Load data
data("e6")

# Create model
model <- create_bvecmodel(e6, p = 4, r = 1,
                          const = "unrestricted", seasonal = "unrestricted",
                          iterations = 10, burnin = 10)
# Number of iterations and burn-in should be much higher.
```
