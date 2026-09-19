# Plot a Spillover Index

Plots the directional connectedness measures of Diebold and Yilmaz
(2012).

## Usage

``` r
# S3 method for class 'bvarspillover'
plot(x, which = "net", ...)
```

## Arguments

- x:

  an object of class 'bvarspillover', usually the result of a call to
  [`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md).

- which:

  which measure to plot. One of `"net"` (default), `"to"` or `"from"`.

- ...:

  further arguments passed to
  [`barplot`](https://rdrr.io/r/graphics/barplot.html).

## Value

`NULL`, invisibly. Called for its side effect.

## Details

One bar per variable, with the credible interval drawn on it where the
object carries one. Net values above zero mark variables that transmit
more forecast error variance than they receive.
