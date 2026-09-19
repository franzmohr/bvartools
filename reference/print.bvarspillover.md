# Print a Spillover Index

Prints the connectedness table of Diebold and Yilmaz (2012).

## Usage

``` r
# S3 method for class 'bvarspillover'
print(x, digits = 1, ...)
```

## Arguments

- x:

  an object of class 'bvarspillover', usually the result of a call to
  [`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md).

- digits:

  the number of decimal places. Defaults to 1, which is what this
  literature reports.

- ...:

  not used.

## Value

`x`, invisibly.

## Details

The body of the table is the posterior mean share of the forecast error
variance of the variable in the row that is attributed to a shock to the
variable in the column, in percent. The `from` column and the `to` row
are the directional spillovers, and the corner is the total index.

Only the means appear here. The credible intervals are in the `total`,
`from`, `to` and `net` elements of the object.
