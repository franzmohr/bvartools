# Prior Inclusion Probabilities

A generic function used to generate prior inclusion probabilities as
required for stochastic search variable selection (SSVS) à la George et
al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
The function invokes particular methods which depend on the class of the
first argument.

## Usage

``` r
inclusion_prior(object, ...)
```

## Arguments

- object:

  an object of a class, for which a method should be called.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`inclusion_prior.bvarmodel`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.bvarmodel.md),
[`inclusion_prior.bvecmodel`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.bvecmodel.md),
[`inclusion_prior.externalforecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md).
