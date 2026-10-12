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

## Examples

``` r

data("e1")
e1 <- diff(log(e1)) * 100
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

summary(model)
#> 
#> Bayesian VAR model with p = 2 
#> 
#> Endogenous variables: invest, income, cons
#> 
#> Variable: invest 
#> 
#>              Mean     SD Naive SD Time-series SD    2.5%     50%    97.5%  
#> invest.l1 -0.2711 0.1078  0.02411        0.01735 -0.4482 -0.2780 -0.08252 *
#> income.l1  0.5708 0.5968  0.13344        0.13344 -0.4470  0.6370  1.63137  
#> cons.l1    0.5750 0.6573  0.14699        0.14699 -0.3868  0.5762  1.68870  
#> invest.l2 -0.1692 0.1193  0.02667        0.01718 -0.3485 -0.1888  0.06431  
#> income.l2  0.1430 0.4339  0.09702        0.09702 -0.5066  0.1102  0.98303  
#> cons.l2    0.6058 0.5336  0.11931        0.11931 -0.4027  0.6886  1.43651  
#> const     -1.1196 1.1307  0.25283        0.25283 -2.8796 -1.1013  0.99097  
#> 
#> Variable: income 
#> 
#>               Mean      SD Naive SD Time-series SD      2.5%      50%   97.5%  
#> invest.l1  0.04005 0.01911 0.004274       0.003178  0.006166  0.03890 0.06916 *
#> income.l1 -0.11710 0.14756 0.032996       0.032996 -0.359045 -0.11064 0.11624  
#> cons.l1    0.31065 0.16515 0.036929       0.036929  0.048693  0.31482 0.54248 *
#> invest.l2  0.06461 0.03001 0.006710       0.004326  0.022215  0.06144 0.11261 *
#> income.l2  0.05554 0.12683 0.028360       0.028360 -0.147869  0.07675 0.28694  
#> cons.l2    0.01852 0.16371 0.036606       0.036606 -0.201894 -0.03244 0.36738  
#> const      1.21179 0.33408 0.074702       0.074702  0.468394  1.22906 1.65075 *
#> 
#> Variable: cons 
#> 
#>                Mean      SD Naive SD Time-series SD       2.5%       50%
#> invest.l1  0.007103 0.01771 0.003960       0.003960 -0.0256691  0.008226
#> income.l1  0.287120 0.13111 0.029317       0.031015  0.1357410  0.265834
#> cons.l1   -0.281874 0.16543 0.036992       0.036992 -0.6192162 -0.255476
#> invest.l2  0.046677 0.02332 0.005215       0.005215  0.0009376  0.042532
#> income.l2  0.388301 0.10644 0.023800       0.023800  0.2181188  0.403598
#> cons.l2   -0.132911 0.15125 0.033821       0.033821 -0.4058702 -0.130951
#> const      1.205422 0.23921 0.053489       0.053489  0.7649076  1.214782
#>              97.5%  
#> invest.l1  0.03528  
#> income.l1  0.56389 *
#> cons.l1   -0.09668 *
#> invest.l2  0.08527 *
#> income.l2  0.57243 *
#> cons.l2    0.16652  
#> const      1.58371 *
#> 
#> Variance-covariance matrix:
#> 
#>                  Mean     SD Naive SD Time-series SD    2.5%     50%  97.5%  
#> invest_invest 20.4489 2.8811  0.64423        0.64423 16.2392 20.7985 24.868 *
#> invest_income  1.0801 0.7720  0.17262        0.17262 -0.2877  0.9602  2.520  
#> invest_cons    1.6747 0.6583  0.14721        0.14721  0.7296  1.6381  2.912 *
#> income_income  1.3906 0.2728  0.06100        0.09226  1.0517  1.3280  1.837 *
#> income_cons    0.7351 0.1843  0.04122        0.04122  0.5053  0.7170  1.151 *
#> cons_cons      1.0970 0.1818  0.04066        0.04066  0.9033  1.0586  1.495 *
#> 
```
