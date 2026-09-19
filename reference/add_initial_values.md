# Add Initial Values of an MCMC Chain

A generic function used to generate initial values of an MCMC chain. The
function invokes particular methods which depend on the class of the
first argument.

## Usage

``` r
add_initial_values(object, ...)
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
[`add_initial_values.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_initial_values.expandingwindow`](https://franzmohr.github.io/bvartools/reference/add_initial_values.expandingwindow.md),
[`add_initial_values.externalforecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`add_initial_values.modellist`](https://franzmohr.github.io/bvartools/reference/add_initial_values.modellist.md).
