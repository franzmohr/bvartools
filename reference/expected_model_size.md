# Expected Size of a Model

A generic function that calculates how much memory a model will take
once its posterior draws are complete, and so roughly how much disk
space it will take when written to a file.

## Usage

``` r
expected_model_size(object, ...)
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

The functions that draw or write a model warn before they start when the
result would exceed the limit set by `options(bvartools.size_warning)`,
in bytes, which is 1e9, one gigabyte, unless set otherwise, and `Inf`
turns the warnings off.
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
checks the size `expected_model_size()` predicts for the draws it is
about to simulate, and
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
checks the size of the object it is about to write. A package that adds
a model class gets both warnings by adding a method to this generic.

## See also

Methods:
[`expected_model_size.bvarmodel`](https://franzmohr.github.io/bvartools/reference/expected_model_size.bvarmodel.md),
[`expected_model_size.bvecmodel`](https://franzmohr.github.io/bvartools/reference/expected_model_size.bvarmodel.md),
[`expected_model_size.expandingwindow`](https://franzmohr.github.io/bvartools/reference/expected_model_size.bvarmodel.md),
[`expected_model_size.modellist`](https://franzmohr.github.io/bvartools/reference/expected_model_size.bvarmodel.md).
