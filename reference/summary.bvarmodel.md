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
#>                 Mean        SD   Naive SD Time-series SD       2.5%        50%
#> invest.l1 -0.2710843 0.1078255 0.02411052     0.01734665 -0.4482211 -0.2780167
#> income.l1  0.5707558 0.5967823 0.13344458     0.13344458 -0.4469869  0.6370360
#> cons.l1    0.5749813 0.6573390 0.14698546     0.14698546 -0.3867653  0.5762362
#> invest.l2 -0.1692296 0.1192756 0.02667084     0.01717632 -0.3485040 -0.1887727
#> income.l2  0.1429599 0.4338994 0.09702285     0.09702285 -0.5065626  0.1101908
#> cons.l2    0.6058185 0.5335870 0.11931367     0.11931367 -0.4026993  0.6886477
#> const     -1.1196186 1.1306912 0.25283024     0.25283024 -2.8796279 -1.1012824
#>                 97.5%  
#> invest.l1 -0.08252382 *
#> income.l1  1.63136914  
#> cons.l1    1.68870280  
#> invest.l2  0.06431349  
#> income.l2  0.98303128  
#> cons.l2    1.43651231  
#> const      0.99097155  
#> 
#> Variable: income 
#> 
#>                  Mean         SD    Naive SD Time-series SD         2.5%
#> invest.l1  0.04004830 0.01911397 0.004274014    0.003177871  0.006166421
#> income.l1 -0.11709568 0.14756062 0.032995558    0.032995558 -0.359045115
#> cons.l1    0.31064535 0.16515143 0.036928983    0.036928983  0.048692618
#> invest.l2  0.06461042 0.03000616 0.006709582    0.004326362  0.022215356
#> income.l2  0.05554184 0.12682837 0.028359687    0.028359687 -0.147869232
#> cons.l2    0.01851658 0.16370656 0.036605901    0.036605901 -0.201894092
#> const      1.21178516 0.33407884 0.074702299    0.074702299  0.468393863
#>                   50%      97.5%  
#> invest.l1  0.03890173 0.06916204 *
#> income.l1 -0.11063978 0.11624027  
#> cons.l1    0.31482153 0.54248012 *
#> invest.l2  0.06143559 0.11260780 *
#> income.l2  0.07675120 0.28694168  
#> cons.l2   -0.03244045 0.36738384  
#> const      1.22906344 1.65074587 *
#> 
#> Variable: cons 
#> 
#>                   Mean         SD    Naive SD Time-series SD          2.5%
#> invest.l1  0.007102718 0.01770754 0.003959527    0.003959527 -0.0256690731
#> income.l1  0.287120007 0.13110965 0.029317008    0.031014990  0.1357410171
#> cons.l1   -0.281874167 0.16543335 0.036992022    0.036992022 -0.6192161774
#> invest.l2  0.046677089 0.02332320 0.005215226    0.005215226  0.0009376257
#> income.l2  0.388300625 0.10643547 0.023799694    0.023799694  0.2181188476
#> cons.l2   -0.132911419 0.15124996 0.033820519    0.033820519 -0.4058702195
#> const      1.205421510 0.23921137 0.053489289    0.053489289  0.7649076336
#>                   50%       97.5%  
#> invest.l1  0.00822566  0.03527898  
#> income.l1  0.26583419  0.56389146 *
#> cons.l1   -0.25547630 -0.09667680 *
#> invest.l2  0.04253207  0.08526897 *
#> income.l2  0.40359828  0.57242932 *
#> cons.l2   -0.13095124  0.16652431  
#> const      1.21478245  1.58370878 *
#> 
#> Variance-covariance matrix:
#> 
#>                     Mean        SD   Naive SD Time-series SD       2.5%
#> invest_invest 20.4488905 2.8810714 0.64422716     0.64422716 16.2392176
#> invest_income  1.0801096 0.7719686 0.17261743     0.17261743 -0.2877492
#> invest_cons    1.6747301 0.6583340 0.14720796     0.14720796  0.7296410
#> income_income  1.3905671 0.2727833 0.06099620     0.09226229  1.0517412
#> income_cons    0.7350975 0.1843345 0.04121845     0.04121845  0.5052571
#> cons_cons      1.0969633 0.1818170 0.04065552     0.04065552  0.9033478
#>                      50%     97.5%  
#> invest_invest 20.7985389 24.867560 *
#> invest_income  0.9602123  2.520291  
#> invest_cons    1.6380683  2.911969 *
#> income_income  1.3280450  1.837461 *
#> income_cons    0.7170029  1.150769 *
#> cons_cons      1.0585583  1.494808 *
#> 
```
