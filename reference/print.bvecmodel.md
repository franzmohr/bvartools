# Printing Model Information

print method for objects of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an object of class 'bvecmodel'.

- digits:

  the number of significant digits to use when printing.

- ...:

  further arguments passed to or from other methods.

## Value

A data frame of the model's specification, as printed, invisibly. It is
the result of
[`get_model_specifications`](https://franzmohr.github.io/bvartools/reference/get_model_specifications.md)
with the columns renamed for display and, where the model has no
unmodelled variables, the columns `m` and `s` dropped.
