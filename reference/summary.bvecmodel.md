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
#>            Mean       SD Naive SD Time-series SD      2.5%       50%   97.5% 
#> l.R    0.002293 0.008388 0.001876       0.001876 -0.009616  0.003415 0.01650 
#> l.Dp  -0.008641 0.028308 0.006330       0.006330 -0.059514 -0.011795 0.02911 
#> const -0.047982 0.071809 0.016057       0.016057 -0.177392 -0.046239 0.06785 
#> 
#> Variable: Dp 
#> 
#>          Mean      SD Naive SD Time-series SD    2.5%     50%   97.5%  
#> l.R    0.3949 0.08424  0.01884        0.05843  0.2429  0.4027  0.5348 *
#> l.Dp  -1.2755 0.08903  0.01991        0.01086 -1.4300 -1.2769 -1.1396 *
#> const -1.8775 0.64621  0.14450        0.47303 -2.8894 -1.8818 -0.7290 *
#> 
#> Variance-covariance matrix:
#> 
#>          Mean      SD Naive SD Time-series SD    2.5%     50%   97.5%  
#> R_R    0.3159 0.03077  0.00688        0.00688  0.2631  0.3196 0.36178 *
#> R_Dp  -0.1273 0.09373  0.02096        0.02096 -0.2688 -0.1369 0.03778  
#> Dp_Dp  3.2362 0.42262  0.09450        0.09450  2.6928  3.1534 4.17869 *
#> 
```
