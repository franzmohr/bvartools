# Add Priors to Bayesian Models

A generic function used to generate prior vectors and matrices. The
function invokes particular methods which depend on the class of the
first argument.

## Usage

``` r
add_priors(object, ...)
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
[`add_priors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`add_priors.expandingwindow`](https://franzmohr.github.io/bvartools/reference/add_priors.expandingwindow.md),
[`add_priors.externalforecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`add_priors.modellist`](https://franzmohr.github.io/bvartools/reference/add_priors.modellist.md).
