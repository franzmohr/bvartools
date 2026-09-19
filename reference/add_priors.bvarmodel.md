# Add Priors to Bayesian Models

Adds prior specifications to a BVAR model, which was produced by
function
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

## Usage

``` r
# S3 method for class 'bvarmodel'
add_priors(object, coef, sigma, varsel = NULL, ...)
```

## Arguments

- object:

  a list of class 'bvarmodel'.

- coef:

  a named list of prior specifications for the coefficients. It has no
  default and must contain at least `v_i` or `minnesota`. Variances are
  specified as precisions, i.e. as inverses of the variances. See
  'Details'.

- sigma:

  a named list of prior specifications for the error term. It has no
  default, and the elements it must contain depend on argument `error`
  of
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).
  See 'Details'.

- varsel:

  a named list of prior specifications for the variable selection
  algorithm. Required if the model was created with `varsel = "ssvs"` or
  `"bvs"`, and not allowed otherwise. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `priors` added, a list with

- `a`:

  the prior of the coefficients: `type`, the \\M \times 1\\ prior means
  `mu` and the \\M \times M\\ prior precision matrix `v_inv`, where
  \\M\\ is the number of coefficients. With variable selection also
  `inprior` and `include`, for SSVS `tau0` and `tau1`, and for TVP
  models `shape` and `rate` of the state variances.

- `psi`:

  for `error = "gamma+covar"` or `"sv+covar"`, the prior of the error
  covariance coefficients with `type`, `mu`, `v_inv` and `varsel`, and
  the variable selection and state variance elements as for `a`.

- `u_sigma`:

  the prior of the error variances with its `type`: `"wishart"` with
  `df` and `scale`, `"gamma"` with `shape` and `rate`, or `"sv"` with
  `mu`, `v_inv`, `shape`, `rate`, `sigma` and `offset`. With a Minnesota
  prior it also holds `u_sigma_inv`, the inverse error covariance matrix
  the prior was scaled with. Not present for `error = "ald"`.

- `u_scale`:

  for `error = "ald"`, the prior of the scales of the asymmetric Laplace
  distribution with `type`, `shape` and `rate`.

## Details

None of the arguments `coef`, `sigma` and `varsel` provides default
hyperparameters: every value that a model needs must be given in the
list it belongs to. A missing required element raises an error, as does
an element that is not listed below.

Argument `coef` can contain the following elements:

- `v_i`:

  a non-negative numeric specifying the prior precision of the
  coefficients, where 0 gives an uninformative prior. Required unless
  `minnesota` is given. It is also required together with `minnesota` if
  `error` is `"gamma+covar"` or `"sv+covar"`, where it sets the prior
  precision of the error covariance coefficients. The precisions of the
  other coefficients are taken from `minnesota` if it is given, and from
  `varsel` for SSVS.

- `v_i_det`:

  a numeric specifying the prior precision of coefficients corresponding
  to deterministic terms. If it is not given, `v_i` is used. Not used if
  `minnesota` is given or SSVS is applied.

- `coint_var`:

  a logical specifying whether the prior mean of the first own lag of an
  endogenous variable should be set to 1, which is commonly used for
  cointegrated VAR models. Default is `FALSE`.

- `const`:

  a numeric or character specifying the prior mean of coefficients,
  which correspond to the intercept. If a numeric is provided, all prior
  means are set to this value. If `coef$const = "mean"`, the mean of the
  respective endogenous variable is used as prior mean. If
  `coef$const = "first"`, the first values of the respective endogenous
  variable is used as prior mean. This can be useful for trend
  estimation when using a TVP model.

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
  `kappa3` will be ignored.

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

Argument `sigma` must contain the elements that belong to the `error` of
the model:

- `"wishart"`: `df` and `scale`. Not available for structural models.

- `"gamma"` and `"gamma+covar"`: `shape` and `rate`.

- `"sv"` and `"sv+covar"`: `mu`, `v_i`, `shape`, `rate`,
  `state_variance` and `offset`.

- `"ald"`: `shape` and `rate`.

The elements are

- `df`:

  a positive integer, or a character expression in `k`, the number of
  endogenous variables, such as `"k"` or `"k + 3"`, specifying the prior
  degrees of freedom of the inverse Wishart prior.

- `scale`:

  a positive numeric specifying the prior error variance of the
  endogenous variables in the inverse Wishart prior.

- `shape`:

  for `"gamma"` and `"gamma+covar"` a non-negative numeric, or a
  character expression in `k` as for `df`, specifying the shape of the
  gamma prior on the error precisions, the inverse error variances,
  whose mean is `shape / rate`. For `"ald"` a positive numeric
  specifying the shape of the inverse gamma prior on the scale of the
  asymmetric Laplace distribution. For both either one value or one per
  endogenous variable. For models with stochastic volatility a numeric
  specifying the shape of the gamma prior on the precision of the state
  equation of the log-volatilities.

- `rate`:

  a positive numeric specifying the rate that corresponds to `shape`,
  for `"gamma"`, `"gamma+covar"` and `"ald"` either one value or one per
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
| \\\frac{\kappa\_{1}}{l}\\ | for own lags of endogenous variables, |
| \\\frac{\kappa\_{2}}{l}\\ | for other endogenous variables, |
| \\\frac{\kappa\_{3}}{1 + l}\\ | for exogenous variables, |
| \\\kappa\_{2}\\ | for contemporaneous endogenous variables of a structural model, |
| \\\kappa\_{4}\\ | for deterministic variables, |

for lag \\l\\ with \\\kappa_1\\, \\\kappa_2\\, \\\kappa_3\\,
\\\kappa_4\\ as the first, second, third and forth element in
`varsel$minnesota`, respectively.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
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
                          iterations = 10, burnin = 10)
# Number of iterations and burn-in should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))
```
