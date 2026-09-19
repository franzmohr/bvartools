# Stochastic Search Variable Selection Prior

Calculates the priors for a Bayesian VAR model, which employs stochastic
search variable selection (SSVS).

## Usage

``` r
# S3 method for class 'bvecmodel'
ssvs_prior(object, tau = c(0.05, 10), semiautomatic = NULL, ...)
```

## Arguments

- object:

  an object of class `"bvarmodel"`, usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- tau:

  a numeric vector of two elements containing the prior standard errors
  of restricted variables (\\\tau_0\\) as its first element and
  unrestricted variables (\\\tau_1\\) as its second. Default is
  `c(0.05, 10)`.

- semiautomatic:

  an optional numeric vector of two elements containing the factors by
  which the standard errors associated with an unconstrained least
  squares estimate of the VAR model are multiplied to obtain the prior
  standard errors of restricted (\\\tau_0\\) and unrestricted
  (\\\tau_1\\) variables. This is the semiautomatic approach described
  in George et al. (2008).

- ...:

  arguments passed forward to method.

## Value

A list containing the vectors of prior standard deviations for
restricted and unrestricted variables, respectively.

## References

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

## Examples

``` r

# Prepare data
data("e6")

# Generate model input
object <- create_bvecmodel(e6, p = 2, r = 1)

# Obtain SSVS prior
prior <- ssvs_prior(object, semiautomatic = c(.1, 10))
```
