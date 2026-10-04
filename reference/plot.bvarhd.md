# Plot a Historical Decomposition

Stacked bars of the contributions of the shocks in each period, and the
part of the response they account for as a line.

## Usage

``` r
# S3 method for class 'bvarhd'
plot(x, baseline = FALSE, ...)
```

## Arguments

- x:

  an object of class `"bvarhd"`, the result of
  [`historical_decomposition`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.md).

- baseline:

  logical: should the baseline be stacked with the shocks? Defaults to
  `FALSE`, which plots the response net of the baseline – the part the
  shocks explain.

- ...:

  further arguments passed to
  [`barplot`](https://rdrr.io/r/graphics/barplot.html).

## Value

`x`, invisibly.
