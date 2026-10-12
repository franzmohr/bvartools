# Single Quantiles of a Quantile Grid

Splits a structural quantile VAR, a model created with
`quantile_grid = TRUE` in
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
into one model per quantile.

## Usage

``` r
split_quantile_grid(object)
```

## Arguments

- object:

  an object of class 'bvarmodel' holding a grid of quantiles, usually
  the result of a call to
  [`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html).

## Value

A list of class 'modellist' with one object of class 'bvarmodel' per
quantile, in increasing order.

## Details

The draws of a grid are stacked level by level, and each level's draws
are the chain the single-quantile model would have drawn. The models
returned hold those draws with the specification of that model, so that
everything written for a single quantile –
[`summary`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md),
[`plot`](https://franzmohr.github.io/bvartools/reference/plot.bvarmodel.md),
[`add_posterior_loglik`](https://rdrr.io/pkg/bayests/man/add_posterior_loglik.html)
– reads them as such. What describes the grid as a whole, its forecasts
and its log likelihood, is not carried over. The impulse responses of a
structural grid are differences of forecasts at a fixed level with and
without a scenario; see `forecast_quantile` in
[`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html).

## See also

[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`add_sign_zero_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`historical_decomposition()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.md),
[`historical_decomposition.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.bvarmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "const",
                           structural = TRUE, error = "ald",
                           quantile = c(0.1, 0.5, 0.9), quantile_grid = TRUE,
                           iterations = 20, burnin = 10)
object <- add_priors(object, coef = list(v_i = 1),
                     sigma = list(shape = 3, rate = .01))
object <- add_initial_values(object)
object <- add_posterior_coefficients(object)

levels <- split_quantile_grid(object)
summary(levels[[3]])
#> 
#> Bayesian Quantile-SVAR model with p = 1 and q = 0.9 
#> 
#> Endogenous variables: Dp, u, r
#> 
#> Period: 194 
#> 
#> Variable: Dp 
#> 
#>          Mean      SD Naive SD Time-series SD      2.5%     50%   97.5%  
#> Dp.l1 0.76345 0.10104 0.022592        0.07446  0.574667 0.79553 0.88183 *
#> u.l1  0.05500 0.02665 0.005959        0.01012  0.002435 0.06057 0.09459 *
#> r.l1  0.06247 0.02273 0.005084        0.01434  0.038808 0.05486 0.10700 *
#> const 0.17492 0.21170 0.047337        0.06311 -0.085387 0.12007 0.59926  
#> Dp    1.00000 0.00000 0.000000        0.00000  1.000000 1.00000 1.00000 *
#> u     0.00000 0.00000 0.000000        0.00000  0.000000 0.00000 0.00000  
#> r     0.00000 0.00000 0.000000        0.00000  0.000000 0.00000 0.00000  
#> 
#> Variable: u 
#> 
#>          Mean      SD Naive SD Time-series SD     2.5%     50%  97.5%  
#> Dp.l1 0.13771 0.10821 0.024196       0.076293 -0.07497 0.16115 0.2735  
#> u.l1  0.99029 0.02146 0.004798       0.009442  0.96136 0.98810 1.0293 *
#> r.l1  0.06062 0.04132 0.009240       0.020163  0.01933 0.04743 0.1515 *
#> const 0.22408 0.11481 0.025671       0.025671  0.08732 0.17866 0.4788 *
#> Dp    0.12034 0.06369 0.014242       0.021003  0.01722 0.10802 0.2307 *
#> u     1.00000 0.00000 0.000000       0.000000  1.00000 1.00000 1.0000 *
#> r     0.00000 0.00000 0.000000       0.000000  0.00000 0.00000 0.0000  
#> 
#> Variable: r 
#> 
#>           Mean      SD Naive SD Time-series SD    2.5%     50%   97.5%  
#> Dp.l1  0.29897 0.15925 0.035610        0.07192  0.0344  0.2876  0.5869 *
#> u.l1   0.79265 0.15687 0.035078        0.08001  0.5466  0.7871  1.0280 *
#> r.l1   0.95891 0.02165 0.004842        0.01293  0.9236  0.9647  0.9904 *
#> const  0.06679 0.15376 0.034382        0.03438 -0.2518  0.1086  0.3032  
#> Dp    -0.76537 0.23629 0.052837        0.14245 -1.2787 -0.7388 -0.4310 *
#> u      0.78200 0.15786 0.035298        0.08916  0.5405  0.7926  1.0141 *
#> r      1.00000 0.00000 0.000000        0.00000  1.0000  1.0000  1.0000 *
#> 
```
