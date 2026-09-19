# Add Initial Values of an MCMC Chain

Adds initial values to a VEC model, which was produced by function
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
in combination with
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

## Usage

``` r
# S3 method for class 'bvecmodel'
add_initial_values(object, method = "maxlik", ...)
```

## Arguments

- object:

  list of class 'bvecmodel'.

- method:

  character specifying the method of how initial values are generated.
  Defaults is `"maxlik"`. Different approaches are used for TVP and SV
  models. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `initial` added, a list with the
starting values `beta` of the cointegration coefficients, a \\K\_\beta r
\times 1\\ matrix, or \\T K\_\beta r \times 1\\ with initial state
`beta_init` for time varying cointegration, and the starting values of
the remaining coefficients and the error term with the elements
described in
[`add_initial_values.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md).

## Details

For argument `method` the following specifications are possible:

- `"maxlik"`:

  Inital values are equal to estimates from maximum likelihood
  regressions.

- `"prior"`:

  Initial values are drawn from the prior. Not possible for
  uninformative priors.

In case `method = "maxlik"`, the initial draw of \\a\\ is the result of
a maximum likelihood estimation of the reduced rank model. When used
with a Wishart prior, the initial draw of \\\Sigma^u\\ is the sum of
squared residuals of ML regression divided by the number of
observations, i.e. \\\frac{uu^{\prime}}{T}\\. In all other cases, the
diagonal elements of \\\Sigma^u\\ are set to the variances of the
variables in \\u\\.

In case `method = "prior"`, all initial draws in the model are random
draws from the respective prior distributions.

In case of a model with time varying parameters (TVP), the initial
states are obtained using the approach specified in argument `method`.
The precisions of the state equations start at the means of their gamma
priors in case `method = "maxlik"` and are drawn from those priors in
case `method = "prior"`. The autocorrelation coefficient `rho` of the
state equation of \\\beta\\ is taken from the prior specification,
whether it is held there or drawn from the prior on it that
`coint$rho_min` and `coint$rho_max` set up.

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
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
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
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 4, r = 1,
                          const = "unrestricted",
                          seasonal = "unrestricted",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)
```
