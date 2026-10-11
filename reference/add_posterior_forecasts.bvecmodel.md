# Add Forecasts

Simulates forecasts of a VEC model, in levels, from its posterior draws.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_posterior_forecasts(
  object,
  forecast_states = NULL,
  scenario = NULL,
  forecast_quantile = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, the result of a call to
  [`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
  and
  [`add_forecast_input`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html).

- forecast_states:

  character, what a VEC model with time varying coefficients or
  stochastic volatility does with them over the forecast horizon.
  `"simulate"` carries each draw's random walks forward, one step per
  period – the loadings and short-run coefficients, the cointegration
  vectors through their state equation, the covariance block and the
  log-volatilities – and rebuilds the VAR in levels from them at every
  period, so that the forecasts are draws from the posterior predictive
  distribution of the estimated model. `"hold"` keeps them at their
  values in the last sample period. If `NULL` (default), the value in
  `object$model$forecast_states` is used, and `"simulate"` when there is
  none. Models with constant coefficients and volatility are unaffected.

- scenario, forecast_quantile:

  not available for a VEC model and refused; see Details.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with `posterior$forecast$forecasts` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and \\Kh\\ columns of forecasts of the levels, stacked by period.
`predict` summarises them. A `forecast_states` that was given is stored
in `model$forecast_states`.

## Details

The forecasts are those of the VAR in levels the VEC model implies. For
a model with constant coefficients they are the forecasts
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
followed by
[`add_posterior_forecasts.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md)
gives. For a model with time varying coefficients they are not: its VAR
representation has no state equation of its own and holds its
coefficients at the last period, while this method steps the states of
the VEC model and converts them to levels anew in every forecast period.

Simulating the cointegration vectors forward depends on the units of the
data. A step \\\eta_t\\ of their state equation moves the error
correction term by \\\eta_t^{\prime} w_t\\, which for series far from
zero – log levels times 100, say – is of the order of their levels
rather than of their variation: the random walk intercept described in
section 'Prior on the cointegration space' of
[`cointspace_prior`](https://franzmohr.github.io/bvartools/reference/cointspace_prior.md),
carried forward with nothing in the data to restrain it. The forecast is
then mostly that drift, and its intervals widen accordingly. A warning
names the series in the error correction term that are more than 50
times their standard deviation per period away from zero at the end of
the sample; `forecast_states = "hold"` is the alternative. Series near
zero, such as interest rates and inflation in percent, do not draw it.

A time varying model estimated on series that
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
centred or scaled cannot be simulated forward and stops with an error
unless `forecast_states = "hold"`: the state equation of its
cointegration vectors belongs to the transformed series, while
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
has put the vectors back on the scale of the data and moved the means
into the constant of the last sample period only.

Simulating the volatility forward needs `posterior$u_sigma_inv$sigma`,
the variance of the log-volatility innovations, which
[`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
stores. A model with stochastic volatility estimated with an earlier
version of the package lacks it and stops with an error unless
`forecast_states = "hold"`.

A model with a seed, see
[`add_seed`](https://rdrr.io/pkg/bayests/man/add_seed.html), draws its
forecast from a seed derived from its own, so two forecasts of it give
the same draws however R's generator stands. A model without one draws
from R's generator as it stands.

A new forecast replaces the forecast errors and predictive scores of the
previous one, which belong to the forecasts they were computed from.

**A VEC model takes no `scenario` and no `forecast_quantile`**: no VEC
sampler conditions its forecast on pinned values, and none estimates a
grid of quantiles. Both are refused rather than ignored, so a model list
mixing VAR and VEC models does not return unconditional forecasts for
its VEC members beside the conditional ones of its VARs.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_seed.bvarmodel.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

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

# Forecasts of the levels, from the VAR in levels the model implies
model <- add_forecast_input(model, n_ahead = 4)
model <- add_posterior_forecasts(model)
predict(model)
#> $fcst
#> , , 1
#> 
#>                R        Dp
#> 1999    3.655367 -1.271014
#> 1999.25 3.085074 -2.750534
#> 1999.5  2.552895 -2.027443
#> 1999.75 2.824815 -3.880204
#> 
#> , , 2
#> 
#>                R        Dp
#> 1999    2.673092 -3.020136
#> 1999.25 1.501671  1.421143
#> 1999.5  1.547385 -2.051648
#> 1999.75 1.882097 -3.590834
#> 
#> , , 3
#> 
#>                R        Dp
#> 1999    4.271828 -6.428048
#> 1999.25 3.890905  3.032415
#> 1999.5  4.098821 -3.864877
#> 1999.75 4.658305  2.148712
#> 
#> , , 4
#> 
#>                R        Dp
#> 1999    4.737414 -3.196932
#> 1999.25 5.007170  1.991408
#> 1999.5  5.324390 -4.171923
#> 1999.75 4.603192  1.583578
#> 
#> , , 5
#> 
#>                R         Dp
#> 1999    3.856632 -0.9303352
#> 1999.25 4.580785  1.8045613
#> 1999.5  4.021135  0.2615851
#> 1999.75 3.638846 -2.2893103
#> 
#> , , 6
#> 
#>                R         Dp
#> 1999    3.731086 -1.3562816
#> 1999.25 3.982091  0.6305008
#> 1999.5  3.300781  0.2365912
#> 1999.75 4.288331 -1.5775657
#> 
#> , , 7
#> 
#>                R          Dp
#> 1999    3.766183 -2.87378762
#> 1999.25 4.060227  0.04331362
#> 1999.5  3.885991  1.25332589
#> 1999.75 4.348217 -1.17391502
#> 
#> , , 8
#> 
#>                R        Dp
#> 1999    3.565818 -3.493270
#> 1999.25 3.053308  3.675179
#> 1999.5  3.612821 -4.520807
#> 1999.75 2.191187  1.147281
#> 
#> , , 9
#> 
#>                R         Dp
#> 1999    3.338786 -0.4907707
#> 1999.25 2.758686  0.7939842
#> 1999.5  2.179595 -1.6797164
#> 1999.75 2.554385  0.4129800
#> 
#> , , 10
#> 
#>                R         Dp
#> 1999    4.193201 -0.5971363
#> 1999.25 4.081350 -0.5478286
#> 1999.5  4.638542 -0.4378599
#> 1999.75 4.896127 -0.4650664
#> 
#> , , 11
#> 
#>                R        Dp
#> 1999    4.386899  1.700942
#> 1999.25 3.577985  1.147479
#> 1999.5  4.024083  1.456880
#> 1999.75 3.574097 -1.195071
#> 
#> , , 12
#> 
#>                R         Dp
#> 1999    3.404396  0.6184863
#> 1999.25 3.422758  1.9191751
#> 1999.5  2.532881 -0.2115203
#> 1999.75 1.796116  0.2389027
#> 
#> , , 13
#> 
#>                R        Dp
#> 1999    3.290271  1.317320
#> 1999.25 3.021228 -2.912794
#> 1999.5  2.583503  2.058786
#> 1999.75 3.512652 -4.382557
#> 
#> , , 14
#> 
#>                R         Dp
#> 1999    4.507624  0.4792462
#> 1999.25 3.488846  0.6529822
#> 1999.5  3.781655  0.9667271
#> 1999.75 2.868592 -1.6141604
#> 
#> , , 15
#> 
#>                R         Dp
#> 1999    3.715280 -0.6014293
#> 1999.25 3.689269 -2.1667627
#> 1999.5  3.364906  4.4015014
#> 1999.75 2.728422 -2.8117527
#> 
#> , , 16
#> 
#>                R         Dp
#> 1999    4.080714  0.4971180
#> 1999.25 4.176871 -0.9897213
#> 1999.5  4.112166 -1.1320953
#> 1999.75 3.738722  2.3682833
#> 
#> , , 17
#> 
#>                R         Dp
#> 1999    3.272646  1.9711302
#> 1999.25 3.360333 -2.3865117
#> 1999.5  4.113068  1.8663872
#> 1999.75 2.908896 -0.9709844
#> 
#> , , 18
#> 
#>                R        Dp
#> 1999    3.924119 -1.825495
#> 1999.25 4.565946  1.923481
#> 1999.5  3.546124  1.926384
#> 1999.75 3.171424 -2.491988
#> 
#> , , 19
#> 
#>                R           Dp
#> 1999    2.982711  2.010761555
#> 1999.25 3.198782 -2.844700052
#> 1999.5  3.439477 -0.002557912
#> 1999.75 3.599131 -2.779284532
#> 
#> , , 20
#> 
#>                R         Dp
#> 1999    3.588307 -0.3848582
#> 1999.25 3.084089  2.5135480
#> 1999.5  3.776709 -0.6006800
#> 1999.75 3.766162  3.6386282
#> 
#> attr(,"tsp")
#> [1] 1999.00 1999.75    4.00
#> 
#> $y
#>            R         Dp
#> 1972 Q3  8.3  1.8871300
#> 1972 Q4  8.7  2.4803600
#> 1973 Q1  8.7  1.6277600
#> 1973 Q2 10.2  0.0289679
#> 1973 Q3  9.8  1.6829000
#> 1973 Q4  9.7  3.8583500
#> 1974 Q1 10.7 -0.1443860
#> 1974 Q2 10.9  1.2723900
#> 1974 Q3 10.8  2.3897200
#> 1974 Q4  9.9  4.2772800
#> 1975 Q1  8.9 -1.0377900
#> 1975 Q2  8.4  0.5584240
#> 1975 Q3  8.7  0.8706090
#> 1975 Q4  8.6  3.8919900
#> 1976 Q1  7.8 -1.9617600
#> 1976 Q2  8.3  0.7319930
#> 1976 Q3  8.1  1.8076400
#> 1976 Q4  7.4  2.3777500
#> 1977 Q1  7.0 -1.4437700
#> 1977 Q2  6.4  1.1005900
#> 1977 Q3  6.0  1.0401700
#> 1977 Q4  6.0  3.5464800
#> 1978 Q1  5.6 -1.5759900
#> 1978 Q2  6.0  1.0019300
#> 1978 Q3  6.4  1.8062100
#> 1978 Q4  6.6  2.5691500
#> 1979 Q1  7.1 -1.7068900
#> 1979 Q2  8.0  0.3504750
#> 1979 Q3  7.8  2.5643800
#> 1979 Q4  8.0  3.1507500
#> 1980 Q1  9.5 -1.3277500
#> 1980 Q2  8.3  1.3096800
#> 1980 Q3  8.3  1.5145800
#> 1980 Q4  9.1  3.1211900
#> 1981 Q1 10.4 -2.0789600
#> 1981 Q2 11.1  1.0448900
#> 1981 Q3 11.3  1.8106000
#> 1981 Q4  9.9  4.0500200
#> 1982 Q1  9.6 -2.0148300
#> 1982 Q2  9.2  0.3808020
#> 1982 Q3  8.8  2.1556900
#> 1982 Q4  8.0  3.2710600
#> 1983 Q1  7.4 -2.0172600
#> 1983 Q2  8.1 -0.3331660
#> 1983 Q3  8.4  1.9287100
#> 1983 Q4  8.3  3.4547300
#> 1984 Q1  7.9 -2.6215600
#> 1984 Q2  8.1 -0.5197530
#> 1984 Q3  7.7  1.3231800
#> 1984 Q4  7.0  3.7427900
#> 1985 Q1  7.7 -2.9194800
#> 1985 Q2  7.0 -0.2613540
#> 1985 Q3  6.4  1.7235300
#> 1985 Q4  6.6  3.7638700
#> 1986 Q1  6.0 -2.0476300
#> 1986 Q2  6.0  0.0900269
#> 1986 Q3  5.8  1.2573200
#> 1986 Q4  6.0  3.4467700
#> 1987 Q1  5.6 -2.2039400
#> 1987 Q2  5.5 -0.2027030
#> 1987 Q3  6.2  0.1734260
#> 1987 Q4  5.8  3.6590100
#> 1988 Q1  5.6 -2.5255700
#> 1988 Q2  6.0  0.0922680
#> 1988 Q3  6.3  0.4513740
#> 1988 Q4  6.2  3.8327200
#> 1989 Q1  7.0 -2.1136800
#> 1989 Q2  7.1 -0.0533104
#> 1989 Q3  7.1  0.9810920
#> 1989 Q4  7.8  3.7354900
#> 1990 Q1  9.0 -1.6232000
#> 1990 Q2  9.0  0.1971240
#> 1990 Q3  9.1  0.6873610
#> 1990 Q4  9.0  2.7458200
#> 1991 Q1  8.6 -1.4078600
#> 1991 Q2  8.6  1.6649700
#> 1991 Q3  8.8  1.4468200
#> 1991 Q4  8.7  4.3004000
#> 1992 Q1  8.2 -1.2776900
#> 1992 Q2  8.4  0.9273530
#> 1992 Q3  8.2  1.7915700
#> 1992 Q4  7.4  3.0508000
#> 1993 Q1  6.5 -1.0590600
#> 1993 Q2  6.7  0.7471080
#> 1993 Q3  6.1  0.4675870
#> 1993 Q4  5.6  3.1067800
#> 1994 Q1  6.2 -1.5368500
#> 1994 Q2  6.9  0.0810146
#> 1994 Q3  7.4  0.8064750
#> 1994 Q4  7.4  2.8765700
#> 1995 Q1  7.1 -1.8378300
#> 1995 Q2  6.4  0.4494670
#> 1995 Q3  6.1  0.9886740
#> 1995 Q4  5.5  2.4594800
#> 1996 Q1  5.8 -1.8518900
#> 1996 Q2  5.9 -0.5901810
#> 1996 Q3  5.5  0.5815030
#> 1996 Q4  5.1  2.4365400
#> 1997 Q1  5.1 -1.5748500
#> 1997 Q2  5.0 -0.7492540
#> 1997 Q3  5.1  0.3882880
#> 1997 Q4  5.1  2.4244800
#> 1998 Q1  4.7 -1.4647000
#> 1998 Q2  4.7 -0.2048970
#> 1998 Q3  4.1  0.2475260
#> 1998 Q4  3.8  2.3923400
#> 
#> attr(,"class")
#> [1] "bvarprd" "list"   
```
