# Stochastic Search Variable Selection Prior

Calculates the priors for a Bayesian VAR model, which employs stochastic
search variable selection (SSVS).

## Usage

``` r
ssvs_prior(object, ...)
```

## Arguments

- object:

  an object of class `"bvarmodel"` or `"bvecmodel"`, usually, a result
  of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- ...:

  arguments passed forward to method.

## Value

A list containing the vectors of prior standard deviations for
restricted and unrestricted variables, respectively.

## References

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

## See also

Methods:
[`ssvs_prior.bvarmodel`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.bvarmodel.md),
[`ssvs_prior.bvecmodel`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.bvecmodel.md),
[`ssvs_prior.externalforecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md).

## Examples

``` r

# Prepare data
data("e1")
data <- diff(log(e1))

# Generate model input
object <- create_bvarmodel(data)

# Obtain SSVS prior
prior <- ssvs_prior(object, semiautomatic = c(.1, 10))
```
