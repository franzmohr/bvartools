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
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

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
[`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md)
– reads them as such. What describes the grid as a whole, its forecasts
and its log likelihood, is not carried over. The impulse responses of a
structural grid are differences of forecasts at a fixed level with and
without a scenario; see `forecast_quantile` in
[`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md).

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
#>             Mean         SD    Naive SD Time-series SD         2.5%        50%
#> Dp.l1 0.69373981 0.05229324 0.011693124    0.021217381  0.594222155 0.69984987
#> u.l1  0.03186382 0.02029465 0.004538021    0.006768237 -0.005073135 0.03358492
#> r.l1  0.05144957 0.01416861 0.003168197    0.008529387  0.028557634 0.05356347
#> const 0.46275990 0.20614535 0.046095501    0.129301563  0.167404127 0.42045853
#> Dp    1.00000000 0.00000000 0.000000000    0.000000000  1.000000000 1.00000000
#> u     0.00000000 0.00000000 0.000000000    0.000000000  0.000000000 0.00000000
#> r     0.00000000 0.00000000 0.000000000    0.000000000  0.000000000 0.00000000
#>            97.5%  
#> Dp.l1 0.77945902 *
#> u.l1  0.06936852  
#> r.l1  0.07425293 *
#> const 0.82260518 *
#> Dp    1.00000000 *
#> u     0.00000000  
#> r     0.00000000  
#> 
#> Variable: u 
#> 
#>             Mean         SD    Naive SD Time-series SD        2.5%        50%
#> Dp.l1 0.17871842 0.10732912 0.023999521    0.023999521 -0.07939719 0.19473590
#> u.l1  0.98762754 0.01765468 0.003947707    0.003947707  0.96530566 0.98486908
#> r.l1  0.03836779 0.01484880 0.003320293    0.004928913  0.02063328 0.03654465
#> const 0.20295331 0.07401143 0.016549460    0.016549460  0.09404142 0.20229349
#> Dp    0.02545400 0.16946292 0.037893060    0.099420548 -0.36369955 0.10282177
#> u     1.00000000 0.00000000 0.000000000    0.000000000  1.00000000 1.00000000
#> r     0.00000000 0.00000000 0.000000000    0.000000000  0.00000000 0.00000000
#>           97.5%  
#> Dp.l1 0.2863116  
#> u.l1  1.0290038 *
#> r.l1  0.0699384 *
#> const 0.3582669 *
#> Dp    0.1865817  
#> u     1.0000000 *
#> r     0.0000000  
#> 
#> Variable: r 
#> 
#>              Mean        SD    Naive SD Time-series SD        2.5%         50%
#> Dp.l1  0.20547267 0.1838521 0.041110580     0.14217832 -0.04844393  0.21581336
#> u.l1   0.59700899 0.1323295 0.029589777     0.02958978  0.39624886  0.58866286
#> r.l1   0.95909850 0.0258854 0.005788152     0.01531496  0.92384267  0.95915443
#> const  0.05037734 0.2811309 0.062862784     0.14309404 -0.27529478 -0.01438054
#> Dp    -0.79745640 0.2695684 0.060277321     0.23653901 -1.17764834 -0.78681710
#> u      0.58910373 0.1435807 0.032105623     0.04853320  0.32700668  0.59353814
#> r      1.00000000 0.0000000 0.000000000     0.00000000  1.00000000  1.00000000
#>            97.5%  
#> Dp.l1  0.4777130  
#> u.l1   0.8912352 *
#> r.l1   0.9966976 *
#> const  0.6651494  
#> Dp    -0.3986366 *
#> u      0.8772428 *
#> r      1.0000000 *
#> 
```
