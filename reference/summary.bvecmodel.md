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

## Examples

``` r

data("e6")
model <- create_bvecmodel(e6 * 100, p = 1, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

summary(model)
#> 
#> Bayesian VEC model with p = 1 
#> 
#> Endogenous variables: R, Dp
#> 
#> Variable: R 
#> 
#>               Mean          SD    Naive SD Time-series SD         2.5%
#> l.R    0.002293456 0.008387707 0.001875548    0.001875548 -0.009615703
#> l.Dp  -0.008641397 0.028307956 0.006329851    0.006329851 -0.059514130
#> const -0.047982398 0.071808753 0.016056925    0.016056925 -0.177391531
#>                50%      97.5% 
#> l.R    0.003414912 0.01649763 
#> l.Dp  -0.011795322 0.02910892 
#> const -0.046239334 0.06785174 
#> 
#> Variable: Dp 
#> 
#>             Mean         SD   Naive SD Time-series SD       2.5%        50%
#> l.R    0.3949499 0.08423987 0.01883661     0.05842869  0.2428555  0.4027062
#> l.Dp  -1.2755401 0.08902769 0.01990720     0.01086221 -1.4299639 -1.2769101
#> const -1.8774965 0.64620867 0.14449665     0.47302626 -2.8894005 -1.8817586
#>            97.5%  
#> l.R    0.5348338 *
#> l.Dp  -1.1395628 *
#> const -0.7290310 *
#> 
#> Variance-covariance matrix:
#> 
#>             Mean         SD    Naive SD Time-series SD       2.5%        50%
#> R_R    0.3159367 0.03076949 0.006880267    0.006880267  0.2630712  0.3196360
#> R_Dp  -0.1273430 0.09373319 0.020959377    0.020959377 -0.2687852 -0.1368659
#> Dp_Dp  3.2361939 0.42262384 0.094501564    0.094501564  2.6927998  3.1533527
#>           97.5%  
#> R_R   0.3617791 *
#> R_Dp  0.0377844  
#> Dp_Dp 4.1786898 *
#> 
```
