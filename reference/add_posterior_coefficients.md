# Posterior Simulation of Model Coefficients

Generic function used for posterior simulation of model coefficients.

## Usage

``` r
add_posterior_coefficients(object, ...)
```

## Arguments

- object:

  an object of a class, for which a method should be called.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## Details

Before the method is called, the size the object will have once its
draws are complete is compared with `options(bvartools.size_warning)`,
and a warning is given if it is larger. See
[`expected_model_size`](https://franzmohr.github.io/bvartools/reference/expected_model_size.md).

## See also

Methods:
[`add_posterior_coefficients.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_coefficients.expandingwindow`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.expandingwindow.md),
[`add_posterior_coefficients.externalforecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`add_posterior_coefficients.modellist`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.modellist.md).
