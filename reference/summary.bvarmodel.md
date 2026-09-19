# Summarising Bayesian VAR Coefficients

summary method for class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
summary(object, ci = 0.95, period = NULL, ...)

# S3 method for class 'summary.bvarmodel'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class 'bvarmodel'.

- ci:

  a numeric between 0 and 1 specifying the probability of the credible
  band. Defaults to 0.95.

- period:

  integer. Index of the period, for which the summary statistics should
  be generated. Only used for TVP or SV models. Default is `NULL`, so
  that the posterior draws of the last time period are used.

- ...:

  further arguments passed to or from other methods.

- x:

  an object of class 'summary.bvarmodel', usually, a result of a call to
  `summary.bvarmodel`.

- digits:

  the number of significant digits to use when printing.

## Value

`summary.bvarmodel` returns a list of class 'summary.bvarmodel', which
contains the following components:

- coefficients:

  A list of various summary statistics of the posterior draws of the VAR
  coefficients.

- sigma:

  A list of various summary statistics of the posterior draws of the
  variance-covariance matrix.

- specifications:

  a list containing information on the model specification.
