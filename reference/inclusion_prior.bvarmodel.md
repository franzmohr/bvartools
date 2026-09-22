# Prior Inclusion Probabilities

Prior inclusion probabilities as required for stochastic search variable
selection (SSVS) à la George et al. (2008) and Bayesian variable
selection (BVS) à la Korobilis (2013).

## Usage

``` r
# S3 method for class 'bvarmodel'
inclusion_prior(
  object,
  prob = 0.5,
  exclude_deterministics = TRUE,
  minnesota_like = FALSE,
  kappa1 = 0.8,
  kappa2 = 0.5,
  kappa3 = 0.5,
  kappa4 = 0.8,
  ...
)
```

## Arguments

- object:

  an object of class `"bvarmodel"`, usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- prob:

  a numeric specifying the prior inclusion probability of all model
  parameters.

- exclude_deterministics:

  logical. If `TRUE` (default), the vector of the positions of included
  variables does not include the positions of deterministic terms.

- minnesota_like:

  logical. If `TRUE`, the prior inclusion probabilities of the
  parameters are calculated in a similar way as the Minnesota prior. See
  'Details'.

- kappa1:

  a numeric specifying the prior inclusion probability of coefficients
  that correspond to own lags of endogenous variables. Only used if
  `minnesota_like = TRUE`. See 'Details'.

- kappa2:

  a numeric specifying the size of the prior inclusion probabilities of
  endogenous variables, which do not correspond to own lags. Only used
  if `minnesota_like = TRUE`. See 'Details'.

- kappa3:

  a numeric specifying the size of the prior inclusion probabilities of
  non-deterministic exogenous variables, between 0 and 1. Default is
  0.5. Only used if `minnesota_like = TRUE`. See 'Details'.

- kappa4:

  a numeric specifying the size of the prior inclusion probabilities of
  deterministic terms. Only used if `minnesota_like = TRUE`. See
  'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

A list containing a matrix of prior inclusion probabilities and an
integer vector specifying the positions of variables, which should be
included in the variable selection algorithm.

## Details

If `minnesota_like = TRUE`, prior inclusion probabilities
\\\underline{\pi}\_1\\ are calculated as

|  |  |
|----|----|
| \\\frac{\kappa_1}{r}\\ | for own lags of endogenous variables, |
| \\\frac{\kappa_2}{r}\\ | for other endogenous variables, |
| \\\frac{\kappa_3}{1 + r}\\ | for unmodelled exogenous variables, |
| \\\kappa_2\\ | for contemporaneous endogenous variables of a structural model, |
| \\\kappa\_{4}\\ | for deterministic variables. |

## Examples

``` r

# Prepare data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model input
object <- create_bvarmodel(e1)

# Obtain inclusion prior
incl <- inclusion_prior(object)
```
