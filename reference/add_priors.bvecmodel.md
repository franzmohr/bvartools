# Add Priors to Bayesian Models

Adds prior specifications to a list of models, which was produced by
function
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

## Usage

``` r
# S3 method for class 'bvecmodel'
add_priors(object, coef, coint, sigma, varsel = NULL, ...)
```

## Arguments

- object:

  a list of class 'bvecmodel'.

- coef:

  a named list of prior specifications for coefficients that do not
  determine the cointegration space. It has no default and must contain
  at least `v_i` or `minnesota`. Variances are specified as precisions,
  i.e. as inverses of the variances. See 'Details'.

- coint:

  a named list of prior specifications for coefficients determining the
  cointegration space of VEC models. It has no default. See 'Details'.

- sigma:

  a named list of prior specifications for the error term. It has no
  default, and the elements it must contain depend on argument `error`
  of
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).
  See 'Details'.

- varsel:

  a named list of prior specifications for the variable selection
  algorithm. Required if the model was created with `varsel = "ssvs"` or
  `"bvs"`, and not allowed otherwise. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `priors` added, a list with

- `beta`:

  the prior of the cointegration space with `type` `"cointspace"`:
  `v_inv` and `p_tau_inv` for constant cointegration parameters, with
  `g_inv` if `coint$g_i` was given, or `rho`, `mu` and `v_inv` of the
  state equation for time varying ones, together with the elements added
  by `p_tau_i = "ml"` or a uniform prior on \\\rho\\.

- `a`:

  the prior of the loadings and the remaining coefficients, with the
  same elements as for a VAR model in
  [`add_priors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md).

- `psi`, `u_sigma`:

  the priors of the error covariance coefficients and error variances,
  as for a VAR model.

## Details

None of the arguments `coef`, `coint`, `sigma` and `varsel` provides
default hyperparameters: every value that a model needs must be given in
the list it belongs to. A missing required element raises an error, as
does an element that is not listed below.

Argument `coef` can contain the following elements:

- `v_i`:

  a non-negative numeric specifying the prior precision of the
  coefficients, where 0 gives an uninformative prior. Required unless
  `minnesota` is given. It is also required together with `minnesota` if
  `error` is `"gamma+covar"` or `"sv+covar"`. The precisions of the
  other coefficients are taken from `minnesota` if it is given, and from
  `varsel` for SSVS.

- `v_i_det`:

  a numeric specifying the prior precision of coefficients corresponding
  to deterministic terms. If it is not given, `v_i` is used. Not used if
  `minnesota` is given or SSVS is applied.

- `const`:

  a numeric or character specifying the prior mean of coefficients,
  which correspond to the intercept. If a numeric is provided, all prior
  means are set to this value. If `coef$const = "mean"`, the mean of the
  respective endogenous variable is used as prior mean. If
  `coef$const = "first"`, the first values of the respective endogenous
  variable is used as prior mean.

- `minnesota`:

  a named list containing the parameters for the calculation of the
  Minnesota prior. It must contain `kappa1`, `kappa2` and `kappa4`, and
  `kappa3` if the model has exogenous variables. For the endogenous
  variable \\i\\ the prior variance of the \\l\\th lag of regressor
  \\j\\ is obtained as \$\$ \frac{\kappa\_{1}}{l^2} \textrm{ for own
  lags of endogenous variables,}\$\$ \$\$ \frac{\kappa\_{1}
  \kappa\_{2}}{l^2} \frac{\sigma\_{i}^2}{\sigma\_{j}^2} \textrm{ for
  endogenous variables other than own lags,}\$\$ \$\$ \frac{\kappa\_{1}
  \kappa\_{3}}{(l+1)^2} \frac{\sigma\_{i}^2}{\sigma\_{j}^2} \textrm{ for
  exogenous variables,}\$\$ \$\$ \kappa\_{1} \kappa\_{4} \sigma\_{i}^2
  \textrm{ for deterministic terms,}\$\$ where \\\sigma\_{i}\\ is the
  residual standard deviation of variable \\i\\ of an unrestricted LS
  estimate. For exogenous variables \\\sigma\_{i}\\ is the sample
  standard deviation. If the model does not contain exogenous variables,
  `kappa3` will be ignored. The function only provides priors for the
  non-cointegration part of the model. However, the residual standard
  errors \\\sigma_i\\ are based on an unrestricted LS regression of the
  endogenous variables on the error correction term and the
  non-cointegration regressors.

- `max_var`:

  a positive numeric specifying the maximum prior variance of the
  coefficients of non-deterministic variables in the Minnesota prior.
  Larger prior variances are set to this value. Only used if `minnesota`
  is given.

- `shape`:

  a numeric specifying the shape of the gamma prior on the precisions,
  the inverse error variances, of the state equation, whose mean is
  `shape / rate`. Required for models with time varying parameters and
  not used otherwise.

- `rate`:

  a numeric specifying the rate of the gamma prior on the precisions of
  the state equation. Required for models with time varying parameters
  and not used otherwise.

- `rate_det`:

  a numeric specifying the prior rate parameter of the error variances
  of the state equation for coefficients, which correspond to
  deterministic terms. If it is not given, `rate` is used. Only used for
  models with time varying parameters.

- `rate_alpha`:

  a numeric specifying the rate of the gamma prior on the precisions of
  the state equation of the loadings. If it is not given, `rate` is
  used. Only used for models with time varying parameters and a positive
  rank. The loadings multiply the levels in the error correction term,
  so a drift that is negligible for a coefficient on a differenced
  regressor moves the fitted value by much more; a rate several orders
  of magnitude below `rate` keeps their drift from absorbing the
  residuals while the other coefficients vary. It does not reach the
  steps of the cointegration vectors, which can absorb the residuals as
  well; see section 'Prior on the cointegration space'.

- `omega_v`:

  a positive numeric, in place of `shape`, `rate`, `rate_det` and
  `rate_alpha`: the variance of a normal prior on the signed standard
  deviation of the state innovations of the coefficients, loadings
  included unless `omega_v_alpha` is given, and of the covariance
  coefficients – the non-centred parameterisation that lets the
  posterior carry the test for time variation of
  [`time_variation_test`](https://franzmohr.github.io/bvartools/reference/time_variation_test.md);
  see
  [`add_priors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md)
  for the model and the draws it adds. The cointegration space keeps its
  state equation. Only for models with time varying parameters,
  `tvp = TRUE`, and `error = "sv"`, `"sv+covar"`, `"gamma"` or
  `"gamma+covar"`.

- `omega_v_alpha`:

  a positive numeric, the `omega_v` of the loadings, as `rate_alpha` is
  their rate under the centred prior. If it is not given, `omega_v` is
  used. The loadings multiply the levels in the error correction term,
  so a standard deviation of their innovations that is small for a
  coefficient on a differenced regressor moves the fitted value by much
  more, and the test for time variation of
  [`time_variation_test`](https://franzmohr.github.io/bvartools/reference/time_variation_test.md)
  compares each loading's posterior at zero with its own prior: a prior
  on the wrong scale makes the Bayes factor of a loading say more about
  the prior than about the loading. Only used with `omega_v` and a
  positive rank.

Argument `coint` specifies the prior on the cointegration space. Its
elements, and the priors they give, are described in section 'Prior on
the cointegration space' below, which is shared with
[`cointspace_prior`](https://franzmohr.github.io/bvartools/reference/cointspace_prior.md).

Argument `sigma` must contain the elements that belong to the `error` of
the model:

- `"wishart"`: `df` and `scale`. Not available for structural models.

- `"gamma"` and `"gamma+covar"`: `shape` and `rate`.

- `"sv"` and `"sv+covar"`: `mu`, `v_i`, `shape`, `rate`,
  `state_variance` and `offset`; with `tvp = TRUE`, `omega_v` may take
  the place of `shape` and `rate`, as `coef$omega_v` does for the
  coefficients.

The elements are

- `df`:

  a positive integer, or a character expression in `k`, the number of
  endogenous variables, such as `"k"` or `"k + 3"`, specifying the prior
  degrees of freedom of the inverse Wishart prior. The samplers add the
  rank \\r\\ of the cointegration matrix to the posterior degrees of
  freedom, as the prior of the loadings requires.

- `scale`:

  a positive numeric specifying the prior error variance of the
  endogenous variables in the inverse Wishart prior.

- `shape`:

  for `"gamma"` and `"gamma+covar"` a non-negative numeric, or a
  character expression in `k` as for `df`, specifying the shape of the
  gamma prior on the error precisions, the inverse error variances,
  whose mean is `shape / rate`, either one value or one per endogenous
  variable. For models with stochastic volatility a numeric specifying
  the shape of the gamma prior on the precision of the state equation of
  the log-volatilities.

- `rate`:

  a positive numeric specifying the rate that corresponds to `shape`,
  for `"gamma"` and `"gamma+covar"` either one value or one per
  endogenous variable.

- `mu`:

  numeric of the prior mean of the initial state of the
  log-volatilities. Only used for models with time varying volatility.

- `v_i`:

  numeric of the prior precision of the initial state of the
  log-volatilities. Only used for models with time varying volatility.

- `state_variance`:

  numeric of the initial draw for the variance of the log-volatilities.
  Only used for models with time varying volatility.

- `offset`:

  numeric of the constant, which is added before taking the log of the
  squared errors. Only used for models with time varying volatility.

For structural models only a gamma prior or stochastic volatility
specification is allowed.

Argument `varsel` can contain the following elements:

- `inprior`:

  a numeric between 0 and 1 specifying the prior probability of a
  variable to be included in the model.

- `covar`:

  logical indicating if the variable selection algorithm should also be
  applied to the error covariance matrix.

- `exclude_det`:

  logical indicating if deterministic terms should be excluded from the
  variable selection algorithm.

- `minnesota`:

  a numeric vector of length 4 containing parameters for the calculation
  of the Minnesota-like inclusion priors. See below.

- `tau`:

  a numeric vector of two elements containing the prior standard errors
  of restricted variables (\\\tau_0\\) as its first element and
  unrestricted variables (\\\tau_1\\) as its second. Only used for SSVS.

- `semiautomatic`:

  an numeric vector of two elements containing the factors by which the
  standard errors associated with an unconstrained least squares
  estimate of the model are multiplied to obtain the prior standard
  errors of restricted (\\\tau_0\\) and unrestricted (\\\tau_1\\)
  variables, respectively. This is the semiautomatic approach described
  in George et al. (2008). Only used for SSVS.

In the case of SSVS, either `tau` or `semiautomatic` must be specified.

If `varsel$minnesota` is specified, prior inclusion probabilities are
calculated in a Minnesota-like fashion as

|  |  |
|----|----|
| \\\frac{\kappa_1}{l}\\ | for own lags of endogenous variables, |
| \\\frac{\kappa_2}{l}\\ | for other endogenous variables, |
| \\\frac{\kappa_3}{1 + l}\\ | for exogenous variables, |
| \\\kappa_2\\ | for contemporaneous endogenous variables of a structural model, |
| \\\kappa\_{4}\\ | for deterministic variables, |

for lag \\l\\ with \\\kappa_1\\, \\\kappa_2\\, \\\kappa_3\\,
\\\kappa_4\\ as the first, second, third and forth element in
`varsel$minnesota`, respectively.

## Prior on the cointegration space

Argument `coint` can contain the following elements. Which of them are
required depends on whether the cointegration vectors are constant or
time varying, since the two are given different kinds of prior: a
matric-variate prior on the cointegration space in the first case, and a
state equation in the second. Any other element raises an error.

- `v_i`:

  non-negative numeric specifying the shrinkage of the cointegration
  space prior, or `"ml"`. See below. Required for models with constant
  cointegration parameters and not used otherwise.

- `p_tau_i`:

  the inverse of the matrix \\P\_\tau\\, which determines the central
  location of the cointegration space \\sp(\beta)\\. Either a numeric of
  its diagonal elements, a full symmetric matrix, or `"ml"`, with one
  row and column per series in the error correction term `data$train$w`.
  See below. Required for models with constant cointegration parameters.
  For models with time varying cointegration parameters only `"ml"` is
  used.

- `weight`:

  positive numeric specifying the weight of a prior centred on the
  maximum likelihood estimate in units of the information of the sample.
  Default is 1. Only used if `p_tau_i = "ml"`.

- `rho`:

  a numeric specifying the autocorrelation coefficient of the state
  equation of \\\beta\\. It must be smaller than 1. Required for models
  with time varying cointegration parameters and not used otherwise. If
  `rho_min` and `rho_max` are given as well, this is the value the chain
  starts \\\rho\\ at rather than the value it keeps, and it must lie
  between them.

- `rho_min`, `rho_max`:

  numerics specifying the support of a uniform prior on \\\rho\\, which
  makes it a drawn parameter rather than a fixed hyperparameter.
  Optional, and either both or neither; they must satisfy \\0 \<
  \\`rho_min`\\ \< \\`rho_max`\\ \le 1\\. Koop et al. (2011) use
  \\(0.999, 1)\\. Only used for models with time varying cointegration
  parameters.

- `g_i`:

  the inverse of the matrix \\G\\ that scales the prior of the loadings,
  for models with constant cointegration parameters and stochastic
  volatility, `error = "sv"` or `"sv+covar"`, and refused for every
  other model. Either a numeric of its diagonal elements, a full
  symmetric positive definite matrix with one row and column per
  endogenous variable, or `"ml"` for the inverse of the maximum
  likelihood estimate of the error covariance. Optional. See below.

For a model with constant cointegration parameters the prior is that of
Koop et al. (2010). The sampler uses `v_i` and `p_tau_i` only through
their product, so with `v_i = 0` the prior on the cointegration space is
uniform whatever `p_tau_i` is. An informative prior on the space
therefore needs a positive `v_i`, which also shrinks the loadings: for
\\\beta\\ close to the centre of the space they have prior \\N(0, \Sigma
/ v)\\.

In Koop et al. (2010) the loadings' prior is scaled by a matrix \\G\\,
which may be the error covariance \\\Sigma\\ or any fixed, known matrix.
The models with a constant error covariance take \\G = \Sigma\\. With
stochastic volatility the covariance differs from period to period, so
\\G\\ is fixed for the whole run instead: \\G^{-1}\\ is `g_i` if it is
given, and otherwise the error precision implied by the starting values
of the log-volatilities, averaged over the sample once before the first
draw. That fallback depends on
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md),
so giving `g_i` makes the prior independent of how the chain is started.
`g_i = "ml"` uses Johansen's (1995) estimate of the error covariance,
the same one `v_i = "ml"` is based on.

With `p_tau_i = "ml"` the prior is centred on the space spanned by
Johansen's (1995) maximum likelihood estimate \\\hat{\beta}\\, computed
from the error correction term as it is stored in the model – so call
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
first if the series should be scaled; scaling is refused afterwards. Let
\\H\\ be an orthonormal basis of that space, \\H\_\perp\\ one of its
orthogonal complement and \\\beta = H + H\_\perp \delta\\. The matrix
\$\$P\_\tau^{-1} = H H^{\prime} + H\_\perp T^{-1} H\_\perp^{\prime},
\quad T = \frac{v}{w} (H\_\perp^{\prime} S\_{11} H\_\perp)^{-1},\$\$
where \\S\_{11}\\ is the cross product of the residuals of a regression
of the error correction term on the short-run regressors and \\w\\ is
`weight`, makes the prior of \\\delta\\ given \\\alpha\\ \\N(0,
(\alpha^{\prime} \Sigma^{-1} \alpha)^{-1} \otimes w^{-1}
(H\_\perp^{\prime} S\_{11} H\_\perp)^{-1})\\. This is the asymptotic
distribution of the maximum likelihood estimator with its precision
multiplied by \\w\\: `weight = 1` adds as much information about the
space as the sample itself holds. Since the centre is estimated from the
same sample, the posterior then overstates the precision of
\\sp(\beta)\\. Eigenvalues of \\T\\ are capped at one, which is the
uniform prior.

With `v_i = "ml"` the shrinkage is chosen so that the prior mean of
\\tr(\alpha^{\prime} \Sigma^{-1} \alpha)\\ is 100 times its maximum
likelihood estimate. It can be combined with a numeric `p_tau_i`.

For a model with time varying cointegration parameters the state
equation is \\\beta_t = \rho \beta\_{t-1} + \eta_t\\ with \\\eta_t \sim
N(0, I)\\, and the prior on the state before the sample is that
equation's own stationary distribution, \\N(0, I / (1 - \rho^2))\\. This
is what makes the prior proper, so a \\\rho\\ close to one is intended
and one further from it draws a warning.
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
gives the loadings the compensating scale: only the product \\\alpha
\beta^{\prime}\\ is identified, so their prior variance is shrunk by
\\1 - \rho^2\\, leaving the product on the scale `coef$v_i` asks for.

When \\\rho\\ is drawn, those last two are computed from `coint$rho`
once and do not follow the draw: the state before the sample keeps the
normal prior built here, and the loadings keep the shrinkage of
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).
This is a deliberate departure from Koop et al. (2011), in whom the
state before the sample is the stationary distribution of whatever
\\\rho\\ currently is. It is also what makes the draw an exact Gibbs
block rather than their Metropolis-within-Gibbs step: with the state
before the sample free of \\\rho\\, the conditional posterior of
\\\rho\\ is a normal truncated to the prior's support. Draws of \\\rho\\
are returned in `object$posterior$beta$rho`.

How far the space moves is not set by the prior on any variance. The
steps \\\eta_t\\ have unit variance whatever `coef$rate` and
`coef$rate_alpha` are, and how far a step moves \\\Pi_t = \alpha_t
\beta_t^{\prime}\\ depends on the scale of \\\beta_t\\. Only the product
is identified, so the posterior can pair small loadings with large
cointegration vectors, which a unit step barely turns, or large loadings
with small ones, which it turns a lot. Two things decide which. One is
the prior precision of the loadings, `coef$v_i` \\/ (1 - \rho^2)\\,
which restrains the drift only if it is tight relative to the scale of
\\\Pi\\ the data call for, and that scale depends on the units of the
data. The other is the series in the error correction term: a step moves
\\\beta_t^{\prime} w_t\\ by \\\eta_t^{\prime} w_t\\, which grows with
the levels of the series and not only with their variation. For log
levels far from zero the shift is of the order of the levels themselves,
and the sampler either lets it act as a random walk intercept, which no
prior on the deterministic terms controls and which can absorb the
residuals of an equation, or draws loadings close to zero, which
switches the error correction term off.

Both have been seen. In the US sub-model of the bgvars data set
`gvar2023`, with the oil price among its endogenous variables and a rank
of one, a prior that pinned the state variances of all coefficients, the
loadings and the constant included, at 1e-14 still left the residual
standard deviation of one equation as low as 0.1 to 0.5 of the maximum
likelihood one in some chains, with which equation and how low depending
on the chain, and a saved chain reported a posterior mean log-likelihood
several hundred above the maximum of the constant coefficient model. The
chains do not move between these outcomes within a few thousand draws,
so a single chain can report any of them. The Austrian sub-model under
the same prior lost fit that the constant coefficient model has, and the
chains that were inspected drew loadings close to zero. On data
simulated from a constant coefficient VEC model, whose posterior under
such a prior should reproduce the maximum likelihood fit, single chains
did so only with the series centred and `coef$v_i = 10`, or around zero
and `coef$v_i = 100`. Around 12, where log levels usually are, a chain
drew loadings close to zero, and around zero with `coef$v_i = 1` one
fitted part of the residuals.

The same mechanism makes
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
risky in a model whose coefficients vary over time. It divides the
series by the standard deviation of their differences, which makes them,
and with them both a drift in a loading and a step of \\\beta_t^{\prime}
w_t\\, larger by the same factor. The coefficient paths can then absorb
the residuals of an equation almost entirely: its residual variance –
or, under stochastic volatility, its volatility path – is reported far
below that of a least squares fit of the same regressors and nearly
constant over the sample. In the model of the vignette on time varying
parameters and stochastic volatility in error correction models,
[`vignette("tvp-sv-vec", package = "bvartools")`](https://franzmohr.github.io/bvartools/articles/tvp-sv-vec.md),
a scaled error correction term did that to one equation in each of four
chains with different seeds, and the unscaled one with a small
`coef$rate` in none.

It is safer to leave the series unscaled, to set `coef$v_i` with the
scale of \\\Pi\\ in mind as well as choosing small rates, and to compare
the residual variances with those of a least squares fit of the same
regressors, in chains with different seeds, for a model that is to be
used. Small rates alone do not guarantee that a model passes that check,
and passing it is what says that neither the coefficients nor the
cointegration vectors fit the residuals.
`scale_error_correction(object, scale = FALSE, centre = TRUE)` centres
the series before posterior simulation, which removes the part of a step
of \\\beta_t^{\prime} w_t\\ that comes from the levels of the series
rather than from their variation, and
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
writes the draws back in terms of the series as they are, so that
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
and the forecasts can use them. It does not reach the other channel, the
scale of the loadings, and it is not a full remedy: in the US sub-model
above, centring alone roughly halved how far the error correction term
and the constant moved over the sample and still left the oil price
equation at about 0.8 of the maximum likelihood residual standard
deviation in two chains. Whether a centred model is right is again what
the comparison with least squares shows.

For a model with time varying cointegration parameters `p_tau_i = "ml"`
centres the marginal prior of the cointegration space on Johansen's
estimate as well, using the informative marginal prior of Koop et al.
(2011, working paper version). The state equation becomes \\\beta_t =
\rho (I_r \otimes P\_\tau) \beta\_{t-1} + \eta_t\\ with \\P\_\tau = H
H^{\prime} + H\_\perp T H\_\perp^{\prime}\\: the part of \\\beta_t\\
along \\sp(H)\\ keeps \\\rho\\, the part off it decays faster, and the
mode of the marginal distribution of \\sp(\beta_t)\\ is \\sp(H)\\ in
every period. \\T\\ is chosen so that the prior spread of the tilt of
\\\beta_t\\ away from \\sp(H)\\ in a period, approximately \\T^\* = (1 -
\rho^2)(I - \rho^2 T^2)^{-1}\\, is the sampling variance of Johansen's
estimator with its precision multiplied by `weight`. The state before
the sample is given the stationary distribution the transition implies.
The transition is stored as element `p_tau`.

\\\rho\\ limits how informative this prior can be: even \\T = 0\\ leaves
the tilt a spread of \\1 - \rho^2\\ per period, and a `weight` asking
for more is floored there with a warning.

A larger `weight` does not always give a tighter posterior. It narrows
the prior spread of the tilt in each period by lowering \\T\\, but \\T\\
is also how much of the tilt carries over from one period to the next.
As \\T\\ approaches zero the tilt of each period becomes independent of
the last, the data of neighbouring periods stop informing it, and the
posterior spread of \\sp(\beta_t)\\ can widen again. On data set E6 with
\\\rho = 0.999\\, for example, `weight = 1`, which gives \\T \approx
0.44\\, produced a tighter posterior than `weight = 100`, which reaches
\\T = 0\\. The useful range of `weight` is the one that keeps \\T\\
clearly above zero, and a \\\rho\\ closer to one widens it.

A `weight` so small that \\T\\ would be the identity in every direction
leaves the noninformative prior unchanged. Both \\P\_\tau\\ and the
prior on the state before the sample are computed at `coint$rho`, and do
not follow the draw if \\\rho\\ is drawn. `coint$v_i` has no counterpart
for these models.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Johansen, S. (1995). *Likelihood-based inference in cointegrated vector
autoregressive models*. Oxford: Oxford University Press.

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224–242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference
in a time varying cointegration model. *Journal of Econometrics,
165*(2), 210–220.
[doi:10.1016/j.jeconom.2011.07.007](https://doi.org/10.1016/j.jeconom.2011.07.007)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

# Load data 
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# The same model with time varying parameters. The cointegration vectors are
# then a state path rather than a draw from a cointegration space prior, so
# 'coint' takes the autocorrelation of that path instead of its shrinkage and
# central location, and 'coef' takes the prior of the state error variances.
model <- create_bvecmodel(e6, p = 4, r = 1, tvp = TRUE,
                          const = "unrestricted",
                          seasonal = "unrestricted",
                          iterations = 10, burnin = 10)

model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10,
                                shape = 3, rate = 0.0001),
                    coint = list(rho = 0.999),
                    sigma = list(df = "k", scale = 1))
```
