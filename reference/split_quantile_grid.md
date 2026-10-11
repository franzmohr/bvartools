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
#>             Mean         SD    Naive SD Time-series SD         2.5%        50%
#> Dp.l1 0.76344771 0.10103540 0.022592203     0.07446193  0.574666716 0.79553085
#> u.l1  0.05499625 0.02664818 0.005958714     0.01011525  0.002434898 0.06057067
#> r.l1  0.06247445 0.02273499 0.005083699     0.01434307  0.038808322 0.05485898
#> const 0.17491908 0.21169635 0.047336742     0.06310642 -0.085386565 0.12007398
#> Dp    1.00000000 0.00000000 0.000000000     0.00000000  1.000000000 1.00000000
#> u     0.00000000 0.00000000 0.000000000     0.00000000  0.000000000 0.00000000
#> r     0.00000000 0.00000000 0.000000000     0.00000000  0.000000000 0.00000000
#>            97.5%  
#> Dp.l1 0.88183164 *
#> u.l1  0.09459216 *
#> r.l1  0.10699599 *
#> const 0.59925897  
#> Dp    1.00000000 *
#> u     0.00000000  
#> r     0.00000000  
#> 
#> Variable: u 
#> 
#>             Mean         SD    Naive SD Time-series SD        2.5%        50%
#> Dp.l1 0.13770570 0.10820689 0.024195795    0.076293444 -0.07497076 0.16115457
#> u.l1  0.99029155 0.02145621 0.004797755    0.009442232  0.96135731 0.98810123
#> r.l1  0.06062479 0.04132126 0.009239714    0.020162668  0.01933366 0.04743353
#> const 0.22407943 0.11480515 0.025671211    0.025671211  0.08731751 0.17866288
#> Dp    0.12033996 0.06369137 0.014241823    0.021002827  0.01721750 0.10802233
#> u     1.00000000 0.00000000 0.000000000    0.000000000  1.00000000 1.00000000
#> r     0.00000000 0.00000000 0.000000000    0.000000000  0.00000000 0.00000000
#>           97.5%  
#> Dp.l1 0.2735186  
#> u.l1  1.0293063 *
#> r.l1  0.1514569 *
#> const 0.4788268 *
#> Dp    0.2307205 *
#> u     1.0000000 *
#> r     0.0000000  
#> 
#> Variable: r 
#> 
#>              Mean         SD    Naive SD Time-series SD        2.5%        50%
#> Dp.l1  0.29896843 0.15925493 0.035610484     0.07191537  0.03439939  0.2875877
#> u.l1   0.79264816 0.15687444 0.035078190     0.08000690  0.54664429  0.7871455
#> r.l1   0.95890750 0.02165361 0.004841895     0.01292860  0.92359571  0.9647326
#> const  0.06679168 0.15376155 0.034382129     0.03438213 -0.25184733  0.1085517
#> Dp    -0.76536639 0.23629223 0.052836548     0.14245008 -1.27869772 -0.7388024
#> u      0.78200164 0.15785842 0.035298215     0.08916475  0.54053089  0.7926273
#> r      1.00000000 0.00000000 0.000000000     0.00000000  1.00000000  1.0000000
#>            97.5%  
#> Dp.l1  0.5869363 *
#> u.l1   1.0280175 *
#> r.l1   0.9903825 *
#> const  0.3032161  
#> Dp    -0.4309915 *
#> u      1.0141046 *
#> r      1.0000000 *
#> 
```
