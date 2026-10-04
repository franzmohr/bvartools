# Test for Time Variation

Generic function used to compute Bayes factors for time variation in the
coefficients and volatilities of a model with time varying parameters.

## Usage

``` r
time_variation_test(object, ...)
```

## Arguments

- object:

  an object with suitable posterior draws passed forward to method.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`time_variation_test.bvarmodel`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md),
[`time_variation_test.bvecmodel`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md),
and for a list of models, the windows of an expanding window comparison
or a folder of stored models,
[`time_variation_test.modellist`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md),
[`time_variation_test.expandingwindow`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md)
and
[`time_variation_test.bvarfolder`](https://franzmohr.github.io/bvartools/reference/folder_steps.md).
