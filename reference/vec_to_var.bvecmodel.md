# Transform a VEC Model to a VAR in Levels

An object of class `'bvecmodel'` is transformed into an object of class
`'bvarmodel'`, which contains the VAR representation of the model in
levels.

## Usage

``` r
# S3 method for class 'bvecmodel'
vec_to_var(object, ...)
```

## Arguments

- object:

  an object of class `'bvecmodel'`, usually, the result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
  optionally already estimated with
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- ...:

  arguments passed forward to method.

## Value

An object of class `'bvarmodel'` with the elements `model` and `data` of
the VAR in levels and, if `object` was estimated, element `posterior`
with the transformed draws, including `loglik` if it was present. Priors
and initial values are not carried over.

## Details

A VEC model and its VAR representation in levels are the same model in
two parameterisations, so a posterior draw of the one is a posterior
draw of the other. The transformation is therefore a change of basis,
which is applied draw by draw: \$\$A_1 = A_0 + \Pi^{y} + \Gamma_1, \quad
A_i = \Gamma_i - \Gamma\_{i - 1}, \quad A_p = -\Gamma\_{p - 1}\$\$ for
the endogenous variables and, analogously, one block later for the
unmodelled, non-deterministic variables \$\$B_0 = \Upsilon_0, \quad B_1
= \Pi^{x} - \Upsilon_0 + \Upsilon_1, \quad B_j = \Upsilon_j -
\Upsilon\_{j - 1}, \quad B_s = -\Upsilon\_{s - 1}.\$\$ \\A_0\\ is the
identity matrix unless the model is structural, in which case it is the
matrix of contemporaneous coefficients. It enters \\A_1\\ because it
multiplies \\\Delta y_t\\ in the VEC model and therefore \\y_t\\ in its
VAR representation, which leaves \\A_0 y\_{t - 1}\\ on the right-hand
side. The contemporaneous coefficients themselves are not affected by
the transformation and are carried over unchanged. Deterministic terms
that entered the model unrestricted carry over unchanged. Those
restricted to the cointegration space become ordinary regressors of the
VAR, entering after the unrestricted ones, and are dropped when the
cointegration rank is zero, in which case they had no effect on the
model. The coefficient transformation itself is performed by the C++
implementation of the package's model library, which is also used to
obtain forecasts of VEC models.

The data matrices are reconstructed from the data of the VEC model, so
the resulting object covers exactly the same periods. Since a VEC model
of lag order \\p\\ and its VAR representation use the same number of
observations, no observations are gained or lost. The levels of the
endogenous variables are recovered from their differences and their
first lag in the error correction term, which requires both to be on the
same scale. If the error correction term was scaled or centred with
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md),
the function therefore expects
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
to have been applied before the transformation.

Note that the covariance matrix of the error term is not affected by the
transformation, since both parameterisations describe the same error
term.

Note also that the elements `priors` and `initial` of the input object
are not carried over, because they describe the coefficients of the VEC
model. Use
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
and
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
if the resulting model should be estimated in its VAR form. The same
applies to the results of variable selection algorithms: a coefficient
of the VAR representation is the sum of multiple coefficients of the VEC
model, so no single inclusion indicator describes it and the draws of
those indicators are dropped.

A model with time varying parameters is a sequence of VEC models, one
per period, and so is its VAR representation: the transformation is
applied to the draws of each period in turn, which yields a VAR model
with time varying parameters whose coefficient draws are the paths of
the coefficients in levels. Its forecasts start from the coefficients of
the last period of the sample, and its impulse responses and variance
decompositions are those of the period given in their argument `period`.
The draws of the variances of the state equations of the VEC
coefficients and of \\\rho\\ are dropped, since they describe how the
coefficients of the VEC model drift and have no counterpart among the
coefficients in levels. For the same reason its forecasts hold the
coefficients and volatilities at the last period rather than simulating
them forward: `model$forecast_states` is set to `"hold"`.

A model with constant coefficients and stochastic volatility has a VAR
representation with the same covariance block and the same
log-volatility random walk, so its forecasts simulate the volatility
forward like those of any other VAR model with stochastic volatility;
see
[`add_posterior_forecasts.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md).
A model of that kind estimated with an earlier version of the package
lacks `posterior$u_sigma_inv$sigma`, the variance of the log-volatility
innovations, and its VAR representation holds the volatility instead.

## References

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md)

## Examples

``` r

# Load data
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# Transform to a VAR model in levels
object <- vec_to_var(model)

# The result can be used like any other VAR model, e.g. for forecasting
object <- add_forecast_input(object, n_ahead = 4)
object <- add_posterior_forecasts(object)
```
