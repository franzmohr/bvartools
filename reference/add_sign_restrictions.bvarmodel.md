# Sign Restrictions

Identifies the shocks of an object of class 'bvarmodel' by the signs of
the impulse responses they produce.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_sign_restrictions(
  object,
  restrictions,
  max_tries = 1000,
  period = NULL,
  ...
)

# S3 method for class 'expandingwindow'
add_sign_restrictions(object, ...)

# S3 method for class 'modellist'
add_sign_restrictions(object, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel', containing posterior draws of the
  coefficients and the error covariance.

- restrictions:

  a data frame of the restrictions, with one row per restriction. See
  'Details'.

- max_tries:

  the number of rotations that are tried per posterior draw before the
  draw is given up on. Defaults to 1000.

- period:

  integer. Index of the period, whose draws should be identified. Only
  used for TVP or SV models. Default is `NULL`, so that the posterior
  draws of the last time period are used.

- ...:

  further arguments passed to or from other methods.

## Value

The object of class 'bvarmodel' with the accepted rotations in element
`q` of its `posterior`, one row per posterior draw, and the
specification of the restrictions in element `sign_restrictions` of its
`model`. Draws for which no admissible rotation was found are recorded
as `NA`.

## Details

The function identifies the structural shocks of the VAR model \$\$y_t =
\sum\_{i=1}^{p} A_i y\_{t - i} + u_t,\$\$ with \\u_t \sim N(0,
\Sigma)\\, by searching for rotations of its Choleski factor whose
impulse responses carry the signs that argument `restrictions` asks for.

Write \\P\\ for the lower triangular Choleski factor of \\\Sigma\\ and
\\Q\\ for an orthogonal matrix. Then \\P Q (P Q)^{\prime} = \Sigma\\ for
every such \\Q\\, so a model rotated by \\Q\\ fits the data exactly as
well as the one the posterior draw describes. The likelihood is
therefore silent about which rotation is the right one, and the
restrictions choose among the models it cannot tell apart. The
identification is in consequence *set valued*: what is obtained is not
one impulse response per posterior draw but the collection of those that
the restrictions admit.

For every posterior draw the function draws rotations uniformly over the
orthogonal group and keeps the first one whose responses satisfy every
restriction. Since a sign restriction cannot distinguish a shock from
its own negative, a column that fails is retried with its sign flipped
before the rotation is discarded. A draw for which no admissible
rotation is found within `max_tries` attempts is dropped from the
identified sample, and [`summary`](https://rdrr.io/r/base/summary.html)
reports how many were. A low acceptance rate is a statement about the
restrictions, not a technicality: it says the model rarely produces the
pattern that was asked of it.

Argument `restrictions` is a data frame with the columns

- `impulse`:

  name of the endogenous variable the shock is named after.

- `response`:

  name of the endogenous variable whose response is restricted.

- `sign`:

  either `1` for a response that must be positive or `-1` for one that
  must be negative.

- `horizon`:

  optional integer of the horizon the restriction applies to, counted
  from zero for the impact period. Defaults to `0`.

Restrictions on several shocks may be combined in one data frame. A
shock that no row mentions is left as it is drawn, since nothing in the
restrictions distinguishes one rotation of it from another; only the
restricted shocks should be interpreted.

Only sign restrictions are supported. Zero restrictions on the impact
responses cannot be imposed by rejection, because the set of rotations
that satisfies them has probability zero, and need the algorithm of
Arias et al. (2018) instead.

The accepted rotations are added to the object as element `q` of its
posterior draws, from where
[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md),
[`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md) and
[`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md)
use them under `type = "sign"`.

The result depends on the state of the random number generator, so
[`set.seed`](https://rdrr.io/r/base/Random.html) is needed to reproduce
it.

## References

Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference
based on structural vector autoregressions identified with sign and zero
restrictions: Theory and applications. *Econometrica, 86*(2), 685-720.
[doi:10.3982/ECTA14468](https://doi.org/10.3982/ECTA14468)

Rubio-Ramirez, J. F., Waggoner, D. F., Zha, T. (2010). Structural vector
autoregressions: Theory of identification and algorithms for inference.
*The Review of Economic Studies, 77*(2), 665-696.
[doi:10.1111/j.1467-937X.2009.00578.x](https://doi.org/10.1111/j.1467-937X.2009.00578.x)

Uhlig, H. (2005). What are the effects of monetary policy on output?
Results from an agnostic identification procedure. *Journal of Monetary
Economics, 52*(2), 381-419.
[doi:10.1016/j.jmoneco.2004.05.007](https://doi.org/10.1016/j.jmoneco.2004.05.007)

## See also

Other post-estimation analysis:
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# A demand shock raises investment and consumption on impact
restrictions <- data.frame(impulse = "invest",
                           response = c("invest", "cons"),
                           sign = c(1, 1),
                           horizon = 0)

set.seed(1234)
model <- add_sign_restrictions(model, restrictions)

# Obtain the identified impulse response
ir <- irf(model, impulse = "invest", response = "cons", type = "sign")
```
