# Time Series Windows

Forwards model input to the same function for individual models.

## Usage

``` r
# S3 method for class 'modellist'
window(x, start = NULL, end = NULL, ...)
```

## Arguments

- x:

  an object of class 'modellist'.

- start:

  the start time of the period of interest.

- end:

  the end time of the period of interest.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'modellist'.
