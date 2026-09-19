# Prior on the Cointegration Space

Checks the specification of the prior on the cointegration space of a
VEC model and builds it, in the form in which
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
stores it as element `priors$beta`.

## Usage

``` r
cointspace_prior(object, coint)
```

## Arguments

- object:

  a list of class 'bvecmodel', or any model with the same layout:
  elements `model$rank` and `model$tvp`, and the regressors
  `data$train$y`, `data$train$w` and `data$train$x`.

- coint:

  a named list of prior specifications for the coefficients determining
  the cointegration space. It has no default. See section 'Prior on the
  cointegration space'.

## Value

`NULL` for a model without cointegration, `model$rank = 0`, whose
`coint` is checked all the same. Otherwise a list with
`type = "cointspace"` and, for constant cointegration parameters,
`v_inv` and `p_tau_inv`, or, for time varying ones, `rho`, `mu` and
`v_inv` of the state equation, together with `rho_min` and `rho_max` for
a uniform prior on \\\rho\\ and the transition `p_tau` added by
`p_tau_i = "ml"`.

## Details

[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
calls this function for VEC models, so it rarely needs to be called
directly. It is exported for packages that build VEC models of their own
layout, such as the sub-models of a global VEC model, whose error
correction term also holds weakly exogenous and global variables.
Calling it gives them the same checks and the same prior.

All dimensions are taken from the regressors, not from the model
specification: the prior has one row and column for each series in the
error correction term `data$train$w`, and one block for each of the
`model$rank` cointegration vectors.

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

For a model with constant cointegration parameters the prior is that of
Koop et al. (2010). The sampler uses `v_i` and `p_tau_i` only through
their product, so with `v_i = 0` the prior on the cointegration space is
uniform whatever `p_tau_i` is. An informative prior on the space
therefore needs a positive `v_i`, which also shrinks the loadings: for
\\\beta\\ close to the centre of the space they have prior \\N(0, \Sigma
/ v)\\.

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

## Examples

``` r

# Load data
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 2, r = 1,
                          const = "unrestricted",
                          iterations = 10, burnin = 10)

# A prior centred on the maximum likelihood estimate of the space
prior <- cointspace_prior(model, coint = list(v_i = 0.01, p_tau_i = "ml"))
prior$p_tau_inv
#>           [,1]     [,2]
#> [1,] 26575.837 5452.663
#> [2,]  5452.663 1119.785
```
