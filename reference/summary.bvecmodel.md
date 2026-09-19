# Summarising Bayesian VEC Coefficients

summary method for class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
summary(object, ci = 0.95, period = NULL, ...)

# S3 method for class 'summary.bvecmodel'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class 'bvecmodel'.

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

  an object of class 'summary.bvecmodel', usually, a result of a call to
  `summary.bvecmodel`.

- digits:

  the number of significant digits to use when printing.

## Value

`summary.bvecmodel` returns a list of class 'summary.bvecmodel', which
contains the following components:

- a:

  A list of various summary statistics of the posterior draws of the VEC
  coefficients.

- sigma:

  A list of various summary statistics of the posterior draws of the
  variance-covariance matrix.

- model:

  a list containing information on the model specification.
