# Printing Model Information

print method for objects of class 'modellist'.

## Usage

``` r
# S3 method for class 'modellist'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an object of class 'modellist'.

- digits:

  the number of significant digits to use when printing.

- ...:

  further arguments passed to or from other methods.

## Value

A data frame with one row per model, as printed, invisibly: the
specifications of
[`get_model_specifications`](https://franzmohr.github.io/bvartools/reference/get_model_specifications.md),
without the columns in which every model agrees.
