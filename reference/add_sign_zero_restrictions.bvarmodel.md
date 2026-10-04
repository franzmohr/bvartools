# Sign and Zero Restrictions

Identifies the shocks of an object of class 'bvarmodel' by the signs of
the impulse responses they produce together with responses that are
restricted to be exactly zero.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_sign_zero_restrictions(
  object,
  restrictions,
  draws = NULL,
  one_sided = FALSE,
  smooth = TRUE,
  period = NULL,
  max_tries = 1,
  ...
)

# S3 method for class 'expandingwindow'
add_sign_zero_restrictions(object, ...)

# S3 method for class 'modellist'
add_sign_zero_restrictions(object, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel', containing posterior draws of the
  coefficients and the error covariance.

- restrictions:

  a data frame of the restrictions, with one row per restriction. See
  'Details'.

- draws:

  integer. The number of draws to resample. Defaults to `NULL`, which
  uses the effective sample size of the importance sampler, the number
  of independent draws it actually produced. Asking for more than that
  is allowed and duplicates draws.

- one_sided:

  logical. Should the numerical derivative behind the importance weights
  be taken on one side only? It is about forty percent faster and
  correspondingly less accurate. Defaults to `FALSE`.

- smooth:

  logical. Should the importance weights be Pareto smoothed? Defaults to
  `TRUE`. `FALSE` uses the raw weights of the paper's Algorithm 3, which
  is what to set to reproduce it exactly.

- period:

  integer. Index of the period, whose draws should be identified. Only
  used for TVP or SV models. Default is `NULL`, so that the posterior
  draws of the last time period are used.

- max_tries:

  integer. The largest number of rotations drawn for each posterior
  draw. Defaults to 1, the algorithm of the paper. More tries find
  rotations where the sign restrictions leave a region too small for a
  single try to hit, and the importance weights account for them. See
  'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object of class 'bvarmodel' with its posterior draws resampled, the
accepted rotations in element `q` of its `posterior`, one row per
resampled draw, and the specification of the restrictions in element
`sign_zero_restrictions` of its `model`, together with `max_tries`, the
number of rotations drawn in all (`tries`), the number of draws that
satisfied the sign restrictions, the effective sample size of the
importance sampler and the share of the total weight held by its single
largest draw, the shape `pareto_k` of the distribution fitted to the
tail of the weights and whether they were `smooth`ed. Element
`sign_restrictions` is set alongside it, so that
[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md),
[`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md) and
[`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md)
read the rotations under `type = "sign"` as they do for
[`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md).

## Details

A sign restriction can be imposed by trying rotations until one of them
carries the signs that were asked for, which is what
[`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md)
does. A zero restriction cannot: the rotations that satisfy one form a
set of probability zero, so no number of tries finds a single member of
it. This function uses the algorithms of Arias, Rubio-Ramirez and
Waggoner (2018) instead, and needs at least one zero restriction to be
worth using – with only signs, `add_sign_restrictions` produces the same
draws more cheaply, and this function says so rather than running.

The restrictions are imposed on \\F(A_0, A\_+)\\, the responses stacked
over the horizons the restrictions mention. That function satisfies
\\F(A_0 Q, A\_+ Q) = F(A_0, A\_+) Q\\ for every orthogonal \\Q\\, so a
zero restriction on the response to the \\j\\th shock is a *linear*
restriction on the \\j\\th column of \\Q\\ once the reduced form is
fixed. The columns are therefore built one at a time, each drawn
uniformly from the sphere in the subspace that the zero restrictions and
the columns already placed leave over. Every rotation the function draws
satisfies the zero restrictions exactly rather than approximately.

That construction does not draw from the posterior conditional on the
zero restrictions but from a distribution that differs from it by a
volume element, which the function computes numerically and divides out.
The result is an importance sample: each draw carries a weight, and the
draws are resampled with replacement according to those weights so that
what the function returns is an ordinary, equally weighted sample. A
draw whose rotation fails the sign restrictions is given weight zero and
so never resampled, and a shock is not retried with its sign flipped,
since the sphere its column is drawn from already covers both of its
signs.

**By default each posterior draw gets one rotation**, as in the paper.
With many sign restrictions that can leave nothing: when one rotation in
tens of thousands satisfies them all, a few thousand posterior draws
produce no admissible rotation, and the function stops. `max_tries`
draws more rotations per draw, but not the way
[`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md)
does. Stopping at the first rotation that satisfies the signs would bias
an importance sampler: the single try is what weights each posterior
draw by the probability \\p\\ that its rotations satisfy the signs, and
a draw whose admissible set is tiny would count as much as one whose set
is large. Instead the rotations are drawn until two of them satisfy the
signs, or until `max_tries` have been drawn. The first of them is kept,
and its weight is multiplied by an unbiased estimate of \\p\\ from the
tries (Girshick, Mosteller and Savage, 1946): \\1 / (N - 1)\\ when the
second came at try \\N\\, and \\1 / T\\ when only one came in all \\T\\
tries. The rotation kept does not depend on how many tries it took, so
the weight stays an unbiased estimate of the paper's and the sample
stays exact. The price is variance: an estimated \\p\\ spreads the
weights, and the effective sample size says by how much. Looking for the
second success also means drawing about twice the rotations a first
success needs.

With more than one try a column whose sign restrictions all hold with
the opposite sign is also flipped rather than rejected. That is exact
too: the proposal does not change when a column changes sign, so the
rotation kept has the same distribution, and the probability that a try
succeeds grows by \\2^s\\, \\s\\ the number of shocks carrying sign
restrictions – the same factor for every draw, which cancels when the
weights are normalised. It cuts the tries needed by that factor, 64 for
six sign restricted shocks. `max_tries = 1` flips nothing: it is the
paper's algorithm and draws exactly what the function drew before the
argument existed.

**The draws that come back are a resample and no longer a chain.** Their
order carries no information, several of them may be copies of the same
original draw, and convergence diagnostics computed on them mean
nothing. How much independent information they carry is the effective
sample size, which [`summary`](https://rdrr.io/r/base/summary.html)
reports.

Argument `restrictions` is a data frame with the columns

- `impulse`:

  name of the endogenous variable the shock is named after.

- `response`:

  name of the endogenous variable whose response is restricted.

- `sign`:

  `1` for a response that must be positive, `-1` for one that must be
  negative, and `0` for one restricted to be exactly zero.

- `horizon`:

  optional horizon the restriction applies to, counted from zero for the
  impact period. `Inf` restricts the long-run response. Defaults to `0`.

**The order of the endogenous variables matters here in a way it does
not for sign restrictions alone.** Columns of \\Q\\ are built in the
order the variables appear, and the \\j\\th column has to keep at least
one dimension after the \\j-1\\ columns before it and its own zero
restrictions have been taken out. A shock carrying many zero
restrictions therefore has to be named early. The function refuses an
ordering that leaves a column with nothing to draw and says which shock
it was.

The result depends on the state of the random number generator, so
[`set.seed`](https://rdrr.io/r/base/Random.html) is needed to reproduce
it. It depends on it in two places rather than one: the rotations are
drawn at random, and so are the matrices that complete each shock's
constraints to a square system. Appendix A.3 of the paper says any draw
of the latter defines a valid algorithm, and that is true – but not that
any draw defines an equally efficient one. Holding a model, its
posterior draws and every rotation fixed, different completions have
moved the effective sample size by a factor of four. A run reporting a
poor effective sample size is therefore worth repeating under a
different seed before the restrictions themselves are blamed.

**The weights are Pareto smoothed before they are used.** Nothing in the
algorithm bounds the ratio of two volume elements, so a draw landing
where the proposal put almost no probability and the target a great deal
carries an enormous weight, the effective sample size collapses, and the
resample is a few copies of that one draw. Pareto smoothing, of Vehtari
et al. (2024), fits a generalised Pareto distribution to the largest few
weights and replaces them by its quantiles. On a model where the raw
weights gave an effective sample size of 4 out of 487 accepted draws,
with one of them holding 46 percent of the total weight, smoothing
returned 138 and the largest share fell to 6 percent; on the same model
under three other seeds, where nothing was wrong, it changed the
effective sample size by less than three percent either way. It buys
that at the cost of a small bias, which is the trade the paper argues
for. `smooth = FALSE` takes the raw weights of Algorithm 3 instead.

The shape \\k\\ of the fitted distribution is worth more than the
smoothing. It estimates how heavy the tail of the weights is, and so
says when the sample cannot be trusted rather than leaving that to be
guessed from an effective sample size: the estimator has a finite
variance only for \\k\\ below one half, and above about 0.7 neither it
nor its effective sample size means much. It is reported by
[`summary`](https://rdrr.io/r/base/summary.html) and kept in element
`pareto_k` of `sign_zero_restrictions`.

**The function warns when the importance sample is not fit to
summarise**: when \\k\\ is 0.7 or above, when the effective sample size
is below twenty, so that the draws returned have no percentiles worth
reading, or when it keeps less than a quarter of the information in the
draws that satisfied the signs. The warning says which repair applies. A
well-behaved tail with a small effective sample size is an efficiency
problem that more candidate draws fix, since the effective sample size
grows roughly in proportion to them. A heavy tail is not fixed that way,
because the trouble is where the proposal put its probability rather
than how much of it was drawn; a different seed, a different variable
ordering or weaker restrictions are the things worth trying.

Applied to a 'modellist' or an 'expandingwindow' the function identifies
each member on its own and returns the collection. Each therefore gets
its own effective sample size, and unless `draws` is given the members
come back with different numbers of draws – which is what it means for
one model, or one window of the sample, to support the restrictions less
well than another. Pass `draws` to give them all the same number.

## References

Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference
based on structural vector autoregressions identified with sign and zero
restrictions: Theory and applications. *Econometrica, 86*(2), 685-720.
[doi:10.3982/ECTA14468](https://doi.org/10.3982/ECTA14468)

Girshick, M. A., Mosteller, F., Savage, L. J. (1946). Unbiased estimates
for certain binomial sampling problems with applications. *The Annals of
Mathematical Statistics, 17*(1), 13-23.

Rubio-Ramirez, J. F., Waggoner, D. F., Zha, T. (2010). Structural vector
autoregressions: Theory of identification and algorithms for inference.
*The Review of Economic Studies, 77*(2), 665-696.
[doi:10.1111/j.1467-937X.2009.00578.x](https://doi.org/10.1111/j.1467-937X.2009.00578.x)

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., Gabry, J. (2024). Pareto
smoothed importance sampling. *Journal of Machine Learning Research,
25*(72), 1-58.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`historical_decomposition()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.md),
[`historical_decomposition.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.bvarmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`split_quantile_grid()`](https://franzmohr.github.io/bvartools/reference/split_quantile_grid.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# Consumption does not respond to the investment shock on impact, and income
# rises. Investment is the first variable of the data set, which is what
# leaves its shock a column to draw: see 'Details' on the ordering.
restrictions <- data.frame(impulse = "invest",
                           response = c("cons", "income"),
                           sign = c(0, 1),
                           horizon = 0)

set.seed(1234)
model <- add_sign_zero_restrictions(model, restrictions)

# Obtain the identified impulse response
ir <- irf(model, impulse = "invest", response = "income", type = "sign")
```
