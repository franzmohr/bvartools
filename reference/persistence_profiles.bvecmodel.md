# Persistence Profiles of a VEC Model

Calculates the persistence profiles of the cointegrating relations of an
object of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
persistence_profiles(object, ...)

# S3 method for class 'bvecmodel'
persistence_profiles(
  object,
  n_ahead = 20,
  ci = 0.95,
  keep_draws = FALSE,
  period = NULL,
  ...
)

# S3 method for class 'modellist'
persistence_profiles(object, ...)

# S3 method for class 'expandingwindow'
persistence_profiles(object, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- ...:

  further arguments passed to or from other methods.

- n_ahead:

  the number of horizons the profile is followed over.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- keep_draws:

  logical. If `FALSE` (default) the draws are summarised by their median
  and the band; if `TRUE` they are returned as they are.

- period:

  the period of a model with time varying coefficients the profile is
  calculated at. Defaults to `NULL`, the last period. Ignored for a
  model whose coefficients are constant.

## Value

A list with one element per cointegrating relation, each a matrix with
one row per horizon and, unless `keep_draws` is `TRUE`, the columns
`median`, `lower` and `upper`.

## Details

The persistence profile of Pesaran and Shin (1996) is the response of a
cointegrating relation to a system-wide shock, scaled to one on impact:
\$\$PP_j(h) = \frac{\beta_j' \Psi_h \Sigma \Psi_h' \beta_j}{\beta_j'
\Sigma \beta_j},\$\$ where \\\Psi_h\\ are the moving average
coefficients of the level VAR the model implies, \\\Psi_0 = I\\, and
\\\Sigma\\ is the covariance of its errors. It starts at one by
construction.

**A profile that does not fall to zero says the relation is not
cointegrating, whatever the rank says.** That is what the statistic is
for: a rank chosen by a test or a criterion is an assertion about how
many stationary combinations exist, and the profile is the check on it.
How quickly the profile falls is the speed of convergence to
equilibrium.

Unlike [`irf`](https://franzmohr.github.io/bvartools/reference/irf.md)
and [`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md),
which a VEC model reaches through
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md),
this statistic cannot be taken from the level VAR alone: `vec_to_var`
drops `beta`, and the profile is a statement about the cointegrating
vectors themselves. The level form is used for the moving average
coefficients and the model's own `beta` for the relations.

**For a model with weakly exogenous variables the profile is partial.**
A cointegrating vector of a VECX model spans the domestic variables and
the weakly exogenous ones, and only the domestic block responds to a
shock to this model. The relation is therefore followed over the part of
it this model governs, which is the whole relation only where the model
has no `exogen`. Where the weakly exogenous variables are endogenous to
a larger system – a global model assembled from sub-models – take the
profile there instead, and the function warns to that effect.

## References

Pesaran, M. H., & Shin, Y. (1996). Cointegration and speed of
convergence to equilibrium. *Journal of Econometrics, 71*(1-2), 117–143.
[doi:10.1016/0304-4076(94)01697-6](https://doi.org/10.1016/0304-4076%2894%2901697-6)

## See also

[`persistence_profiles`](https://franzmohr.github.io/bvartools/reference/persistence_profiles.md)
for the generic.

## Examples

``` r

# Load data
data("e6")

# Create model
model <- create_bvecmodel(e6, p = 2, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 0.0001))

model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

profiles <- persistence_profiles(model, n_ahead = 12)
round(profiles[[1]][1:5, ], 3)
#>   median lower upper
#> 0  1.000 1.000 1.000
#> 1  0.170 0.099 0.357
#> 2  0.112 0.001 0.176
#> 3  0.075 0.027 0.219
#> 4  0.003 0.000 0.020
```
