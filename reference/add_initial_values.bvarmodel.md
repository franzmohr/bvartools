# Add Initial Values of an MCMC Chain

Adds initial values to a VAR model.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_initial_values(object, method = "ols", ...)
```

## Arguments

- object:

  list of class 'bvarmodel', usually, the result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  in combination with
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

- method:

  character specifying the method of how initial values are generated.
  Default is `"ols"`. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `initial` added, a list of
starting values with

- `a`:

  the coefficients, an \\M \times 1\\ matrix, or \\TM \times 1\\ for TVP
  models, whose initial states and state precisions are held in `a_init`
  and `a_sigma_inv`.

- `u_sigma_inv`:

  for `error = "wishart"`, the \\K \times K\\ inverse error covariance
  matrix.

- `u_omega_inv`:

  for gamma priors, the \\K \times K\\ diagonal matrix of error
  precisions.

- `h`, `h_init`:

  for stochastic volatility, the \\T \times K\\ log-volatilities and
  their initial states.

- `psi`:

  for `"gamma+covar"` and `"sv+covar"`, the error covariance
  coefficients.

- `a_lambda`:

  with variable selection, the inclusion indicators.

- `u_scale`, `w`:

  for `error = "ald"`, the \\K \times 1\\ scales and the \\T \times K\\
  latent weights of the asymmetric Laplace distribution.

Elements that do not apply to a model are absent.

## Details

For argument `method` the following specifications are possible:

- `"ols"`:

  Inital values are equal to estimates from an unrestricted LS
  regression.

- `"prior"`:

  Initial values are drawn from the prior. Not possible for
  uninformative priors.

In case `method = "ols"`, the initial draw of \\a\\ is the result of the
unrestricted LS regression \\(Z^{\prime}Z)^{-1}Z^{\prime}y\\ of the
model \\y\_{t} = Z\_{t} a + u\_{t}\\, where \\Z\\ is the matrix of
regressors in SUR form and \\y\\ is the vector of endogenous variables.

The initial draw of \\\Sigma^u\\ is the sum of squared residuals divided
by the number of observations, i.e. \\\frac{uu^{\prime}}{T}\\, with
\\u\\ as the residuals of the LS regression.

In case of a model with time varying parameters (TVP), the initial
states are obtained using the approach specified in argument `method`.
The precisions of the state equations start at the means of their gamma
priors in case `method = "ols"` and are drawn from those priors in case
`method = "prior"`. So `method = "ols"` uses no random numbers, and a
seed set after `add_initial_values` fixes the posterior draws.

## Seed

The function also stores the seed of the posterior simulation as element
`seed` of `object$model`, unless the model has one already. It is drawn
from R's random number generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before this call
makes it reproducible.
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
replaces it.

## See also

Other model set-up:
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
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
                          iterations = 50, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)
```
