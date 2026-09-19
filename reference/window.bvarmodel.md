# Time Series Windows

Restricts the observations in the training sample of a model of class
'bvarmodel' to the specified time window.

## Usage

``` r
# S3 method for class 'bvarmodel'
window(x, start = NULL, end = NULL, ...)
```

## Arguments

- x:

  an object of class 'bvarmodel'.

- start:

  the start time of the period of interest.

- end:

  the end time of the period of interest.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'bvarmodel'.

## Details

Posterior draws that form a path over the periods of the training sample
are cut to the periods that remain, so that they still refer to the same
observations as the data. These are the coefficients of a model with
time varying parameters, the draws of the error term that vary by period
– as under stochastic volatility – and the pointwise log-likelihood.
Draws that do not vary by period are left as they are.
