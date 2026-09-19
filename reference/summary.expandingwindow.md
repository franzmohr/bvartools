# Summarising Bayesian Models

summary method for class 'expandingwindow'.

## Usage

``` r
# S3 method for class 'expandingwindow'
summary(object, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow'.

- ...:

  further arguments passed to or from other methods.

## Value

A list of class 'summary.bvarmodel' or 'summary.bvecmodel'.

## Details

The method passes the last element of the expanding window model list to
its corresponding summary function; assuming that this object contains
the model with the highest amount of available observations.
