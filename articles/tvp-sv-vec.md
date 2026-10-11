# Time Varying Parameters and Stochastic Volatility in Error Correction Models

## Introduction

The companion vignette on TVP-SV-VAR models describes what the arguments
`tvp` and `error` do and which priors they require. They mean the same
thing for a vector error correction model, and `create_bvecmodel` takes
them in the same way. What changes is the object they are applied to,
and that change is not cosmetic: in a VEC model the time varying block
includes the cointegration term, so the model allows the long-run
relationship itself to move.

The estimated model is

``` math
\Delta y_t = \Pi_t w_t + \sum_{l=1}^{p-1} \Gamma_{l,t} \Delta y_{t-l} + C_t d_t + u_t,
\qquad u_t \sim N(0, \Sigma_t),
```

with $`\Pi_t = \alpha_t \beta_t^{\prime}`$ of rank $`r`$. The
coefficients of the non-cointegration part and the loadings $`\alpha_t`$
follow random walks and $`\Sigma_t`$ is the stochastic volatility
specification of the VAR vignette. The cointegration vectors get a state
equation of their own,

``` math
\beta_t = \rho \beta_{t-1} + \eta_t, \qquad \eta_t \sim N(0, I),
```

which is the specification of Koop, León-González and Strachan (2011). A
constant coefficient VEC model instead puts a matric-variate prior on
the cointegration space, so `add_priors` asks for different arguments in
the two cases: `coint$v_i` and `coint$p_tau_i` for a constant model,
`coint$rho` for a time varying one.

## Data

Data set `at_macrodata` contains, among the domestic series of Austria,
a short-term interest rate (`r`), a long-term interest rate (`lr`) and
inflation (`Dp`). The model uses the three series in percent per quarter
from 1979Q2 to 2019Q4 and asks the question this model is for, and which
Koop, León-González and Strachan (2011) ask of UK data: whether the
long-run relationship between the two rates and inflation is the same at
the end of the sample as it was at its beginning.

``` r

library(bvartools)

data("at_macrodata")

data <- window(at_macrodata[["domestic"]][, c("r", "lr", "Dp")], end = c(2019, 4)) * 100

plot(data, main = "Austrian interest and inflation rates")
```

![plot of chunk data](figures/tvp-sv-vec-data-1.png)

plot of chunk data

## Model set-up

``` r

model <- create_bvecmodel(data, p = 2, r = 1, const = "unrestricted",
                          tvp = TRUE, error = "sv+covar", varsel = "bvs",
                          iterations = 3000, burnin = 1000)

model[["model"]][["algorithm"]]
#> [1] "VecTvpStochvol"
```

`p` is the lag order of the level VAR, so the model contains
$`p - 1 = 1`$ lag of the differenced series. The cointegration rank is
fixed at one, which for a term structure is the natural specification.

Models with a prior on the cointegration space are sensitive to the
scale of the series in the error correction term, and
`scale_error_correction` can divide them by the standard deviation of
their differences. This model is estimated on the series as they are,
for a reason explained with the priors below, so nothing has to be
rescaled after estimation either.

## Priors

``` r

model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10,
                                shape = 3, rate = 0.000001, rate_det = 0.01),
                    coint = list(rho = 0.999),
                    sigma = list(shape = 3, rate = 0.01,
                                 mu = 0, v_i = 0.01,
                                 state_variance = 0.05, offset = 1e-8),
                    varsel = list(inprior = 0.5, exclude_det = TRUE))
```

Arguments `coef` and `sigma` are the ones of the VAR vignette: prior
precisions for the initial state of the coefficients, an inverse gamma
prior for the variances of their state equation, and the prior of the
stochastic volatility block.

`coef$rate` is two orders of magnitude smaller than in the VAR vignette,
and the series in the error correction term are left unscaled. Both
choices address the same problem. The coefficient paths can absorb most
of the residuals of the short-term interest rate equation, whose
residual standard deviation is then reported at a small fraction of the
one of a least squares fit of the same regressors and almost flat over
forty years – the symptom the VAR vignette describes. In four chains of
the length used here, each with a different seed, that happened in every
chain with `rate = 0.0001`, which put the residual standard deviation of
that equation at between 4 and 53 percent of the least squares value,
and in every chain with `rate = 0.000001` when the series were scaled,
at between 3 and 10 percent. With `rate = 0.000001` and unscaled series,
all four chains put it at about 95 percent.

Neither choice reaches all of the drift, though. The cointegration
vectors take steps of unit variance whatever `coef$rate` is, and how far
those steps move $`\Pi_t`$ depends on the prior precision of the
loadings relative to the scale of $`\Pi`$ and on the levels of the
series in the error correction term, as section ‘Prior on the
cointegration space’ of
[`?add_priors.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md)
explains. Small rates alone therefore do not guarantee that a model
passes the comparison with least squares, which is why that comparison
is the check to make, whichever block would fit the residuals.

Argument `coint` is where a VEC model differs. `rho` is the
autocorrelation of the state equation of $`\beta`$ and must be smaller
than one. The prior on the state before the sample is that equation’s
own stationary distribution, $`N(0, I / (1 - \rho^2))`$, which is what
makes the prior proper, so a value close to one is intended and a value
far from it draws a warning. Only the product
$`\alpha_t \beta_t^{\prime}`$ is identified, and the prior variance of
the loadings is shrunk by $`1 - \rho^2`$ to compensate, which leaves the
product on the scale that `coef$v_i` asks for.

The loadings are never subject to variable selection, so the
cointegration term is always in the model. Variable selection applies to
the coefficients of the differenced regressors, and `exclude_det = TRUE`
keeps the unrestricted constant out of it as well.

## Initial values and posterior draws

Initial values of the state variances and volatilities are drawn from
their priors, so the seed is set before them.

``` r

set.seed(20260912)
model <- add_initial_values(model)
```

``` r

model <- add_posterior_coefficients(model)
```

## Evaluation

### Summary statistics for a period

As for a TVP-SV-VAR, a summary table is the summary of a period, and the
last one is used by default. The first block of columns is the
cointegration matrix $`\Pi_t`$ rather than $`\alpha_t`$ and $`\beta_t`$
separately, since the factors are identified only up to an invertible
$`r \times r`$ transformation while their product is not.

``` r

summary(model)
#> 
#> Bayesian TVP-SV-VEC model with p = 2 
#> 
#> Variable selection algorithm: Bayesian variable selection (Korobilis, 2013)
#> 
#> Endogenous variables: r, lr, Dp
#> 
#> Period: 161 
#> 
#> Variable: r 
#> 
#>                  Mean          SD     Naive SD Time-series SD        2.5%
#> l.r      -0.001848429 0.019240098 0.0003512745   0.0003992379 -0.04400914
#> l.lr     -0.002341233 0.020134146 0.0003675975   0.0004830594 -0.04946484
#> l.Dp      0.007965832 0.046001288 0.0008398648   0.0014310419 -0.08317439
#> d.r.l01   0.002956491 0.030085886 0.0005492906   0.0010748096 -0.01888376
#> d.lr.l01  0.202364073 0.076911049 0.0014041972   0.0043759480  0.00000000
#> d.Dp.l01 -0.001789576 0.008846161 0.0001615081   0.0008554777 -0.03693315
#> const     0.013097995 0.023981803 0.0004378458   0.0009208237 -0.03839760
#>                    50%      97.5% Incl. prob. 
#> l.r      -0.0003989906 0.03718636  1.00000000 
#> l.lr     -0.0004764675 0.04044824  1.00000000 
#> l.Dp      0.0068156107 0.10483619  1.00000000 
#> d.r.l01   0.0000000000 0.08725826  0.08233333 
#> d.lr.l01  0.2067874132 0.34234393  0.96533333 
#> d.Dp.l01  0.0000000000 0.00000000  0.05666667 
#> const     0.0146934733 0.05779979  1.00000000 
#> 
#> Variable: lr 
#> 
#>                   Mean          SD     Naive SD Time-series SD        2.5%
#> l.r      -1.048360e-03 0.025729465 4.697536e-04   5.621982e-04 -0.06045351
#> l.lr     -7.936845e-03 0.029745478 5.430756e-04   1.435499e-03 -0.08242263
#> l.Dp      1.677572e-02 0.063645912 1.162010e-03   1.738439e-03 -0.11298987
#> d.r.l01   2.117568e-03 0.016714946 3.051718e-04   6.229977e-04  0.00000000
#> d.lr.l01  2.796428e-02 0.071553362 1.306380e-03   5.072217e-03  0.00000000
#> d.Dp.l01 -6.403832e-05 0.002118884 3.868536e-05   5.198422e-05  0.00000000
#> const    -2.124139e-02 0.043697645 7.978062e-04   1.310558e-03 -0.10771068
#>                    50%      97.5% Incl. prob. 
#> l.r      -0.0000500699 0.05239768 1.000000000 
#> l.lr     -0.0020780602 0.04031086 1.000000000 
#> l.Dp      0.0148151450 0.14794691 1.000000000 
#> d.r.l01   0.0000000000 0.04725873 0.075333333 
#> d.lr.l01  0.0000000000 0.25816385 0.197666667 
#> d.Dp.l01  0.0000000000 0.00000000 0.008666667 
#> const    -0.0203602721 0.06214006 1.000000000 
#> 
#> Variable: Dp 
#> 
#>                 Mean         SD     Naive SD Time-series SD        2.5%
#> l.r       0.02866077 0.34263081 0.0062555542    0.020175861 -0.66064268
#> l.lr      0.19784602 0.30161969 0.0055067970    0.025502965 -0.38650661
#> l.Dp     -0.83470370 0.23211174 0.0042377613    0.011922120 -1.31360281
#> d.r.l01   0.47786606 0.32720551 0.0059739279    0.016750042  0.00000000
#> d.lr.l01  0.01663549 0.17614079 0.0032158761    0.004014099 -0.38438498
#> d.Dp.l01  0.01288475 0.04675904 0.0008536994    0.002532061  0.00000000
#> const     0.34086449 0.13958206 0.0025484080    0.007434376  0.06485793
#>                  50%      97.5% Incl. prob.  
#> l.r       0.04382958  0.6659260   1.0000000  
#> l.lr      0.19752692  0.8024446   1.0000000  
#> l.Dp     -0.82270787 -0.4112662   1.0000000 *
#> d.r.l01   0.50983476  1.0718300   0.8060000  
#> d.lr.l01  0.00000000  0.4850410   0.2723333  
#> d.Dp.l01  0.00000000  0.1769038   0.1433333  
#> const     0.33906181  0.6103480   1.0000000 *
#> 
#> Variance-covariance matrix:
#> 
#>               Mean           SD     Naive SD Time-series SD         2.5%
#> r_r   0.0004424237 0.0009411751 1.718343e-05   4.295791e-05 6.110585e-07
#> r_lr  0.0001009317 0.0001920165 3.505725e-06   1.016780e-05 1.274708e-07
#> r_Dp  0.0005157273 0.0010637404 1.942115e-05   5.017812e-05 6.073661e-07
#> lr_lr 0.0027242863 0.0011591587 2.116325e-05   9.604966e-05 9.962846e-04
#> lr_Dp 0.0029118933 0.0018269186 3.335482e-05   9.827043e-05 1.852222e-04
#> Dp_Dp 0.0513290548 0.0185101136 3.379469e-04   1.017208e-03 2.411791e-02
#>                50%        97.5%  
#> r_r   2.231408e-04 0.0020215943 *
#> r_lr  4.853952e-05 0.0004861096 *
#> r_Dp  2.503330e-04 0.0023855942 *
#> lr_lr 2.576446e-03 0.0054485960 *
#> lr_Dp 2.615767e-03 0.0070700934 *
#> Dp_Dp 4.847440e-02 0.0951584246 *
```

``` r

period_of <- function(object, year, quarter) {
  which(stats::time(object[["data"]][["train"]][["y"]]) == year + (quarter - 1) / 4)
}

summary(model, period = period_of(model, 1985, 1))
#> 
#> Bayesian TVP-SV-VEC model with p = 2 
#> 
#> Variable selection algorithm: Bayesian variable selection (Korobilis, 2013)
#> 
#> Endogenous variables: r, lr, Dp
#> 
#> Period: 22 
#> 
#> Variable: r 
#> 
#>                  Mean          SD     Naive SD Time-series SD         2.5%
#> l.r      -0.024510001 0.029120490 0.0005316650    0.002981496 -0.096477722
#> l.lr     -0.029279724 0.029492119 0.0005384500    0.002977153 -0.103557632
#> l.Dp      0.113811013 0.061782139 0.0011279824    0.002135183 -0.004177071
#> d.r.l01   0.002944663 0.029981331 0.0005473817    0.001076070 -0.019801107
#> d.lr.l01  0.202371107 0.077048369 0.0014067043    0.004379585  0.000000000
#> d.Dp.l01 -0.001814200 0.008923727 0.0001629242    0.000868795 -0.037130585
#> const    -0.022459142 0.065090082 0.0011883769    0.005643279 -0.152284444
#>                  50%      97.5% Incl. prob. 
#> l.r      -0.01845502 0.01961815  1.00000000 
#> l.lr     -0.02409280 0.01415477  1.00000000 
#> l.Dp      0.11275107 0.23919694  1.00000000 
#> d.r.l01   0.00000000 0.08297552  0.08233333 
#> d.lr.l01  0.20731106 0.34322304  0.96533333 
#> d.Dp.l01  0.00000000 0.00000000  0.05666667 
#> const    -0.02261960 0.10686633  1.00000000 
#> 
#> Variable: lr 
#> 
#>                   Mean          SD     Naive SD Time-series SD        2.5%
#> l.r      -7.100208e-03 0.012467029 2.276158e-04   7.238598e-04 -0.03629710
#> l.lr     -1.008772e-02 0.014564890 2.659173e-04   1.239595e-03 -0.04498621
#> l.Dp      3.722972e-02 0.037213669 6.794255e-04   7.948822e-04 -0.03443994
#> d.r.l01   2.197954e-03 0.016585911 3.028159e-04   6.662399e-04  0.00000000
#> d.lr.l01  2.803078e-02 0.071644324 1.308040e-03   5.091790e-03  0.00000000
#> d.Dp.l01 -7.776189e-05 0.002219762 4.052712e-05   5.075608e-05  0.00000000
#> const    -9.455487e-03 0.036147879 6.599670e-04   1.348269e-03 -0.08080780
#>                   50%      97.5% Incl. prob. 
#> l.r      -0.004736040 0.01396778 1.000000000 
#> l.lr     -0.006851734 0.01107162 1.000000000 
#> l.Dp      0.037153433 0.11059108 1.000000000 
#> d.r.l01   0.000000000 0.05100137 0.075333333 
#> d.lr.l01  0.000000000 0.25817701 0.197666667 
#> d.Dp.l01  0.000000000 0.00000000 0.008666667 
#> const    -0.009573793 0.06266561 1.000000000 
#> 
#> Variable: Dp 
#> 
#>                 Mean        SD     Naive SD Time-series SD       2.5%
#> l.r       0.20163356 0.1932672 0.0035285598    0.022481920 -0.1833285
#> l.lr      0.26913654 0.2088513 0.0038130852    0.034326130 -0.1027361
#> l.Dp     -1.02154926 0.1769649 0.0032309230    0.008685228 -1.4089087
#> d.r.l01   0.47757712 0.3271194 0.0059723567    0.016742929  0.0000000
#> d.lr.l01  0.01663580 0.1761289 0.0032156588    0.004013375 -0.3838768
#> d.Dp.l01  0.01289159 0.0466729 0.0008521266    0.002348725  0.0000000
#> const     0.02008127 0.2871357 0.0052423575    0.047240068 -0.5602479
#>                  50%      97.5% Incl. prob.  
#> l.r       0.20816605  0.5624278   1.0000000  
#> l.lr      0.25804082  0.7002838   1.0000000  
#> l.Dp     -1.01415607 -0.6951284   1.0000000 *
#> d.r.l01   0.51133094  1.0712156   0.8060000  
#> d.lr.l01  0.00000000  0.4874602   0.2723333  
#> d.Dp.l01  0.00000000  0.1761074   0.1433333  
#> const     0.03632376  0.5700426   1.0000000  
#> 
#> Variance-covariance matrix:
#> 
#>              Mean          SD     Naive SD Time-series SD        2.5%
#> r_r   0.027121225 0.016836170 3.073850e-04   0.0006385262 0.007721146
#> r_lr  0.006149327 0.004163403 7.601299e-05   0.0001511844 0.001388687
#> r_Dp  0.031995143 0.022415127 4.092424e-04   0.0008642580 0.007055610
#> lr_lr 0.003874753 0.001507750 2.752762e-05   0.0001235440 0.001788588
#> lr_Dp 0.009767837 0.005775801 1.054512e-04   0.0002130866 0.003048780
#> Dp_Dp 0.164504174 0.051028038 9.316403e-04   0.0033452324 0.090111642
#>               50%       97.5%  
#> r_r   0.023037869 0.072000643 *
#> r_lr  0.005155673 0.016924510 *
#> r_Dp  0.026573953 0.088195955 *
#> lr_lr 0.003567705 0.007465279 *
#> lr_Dp 0.008326149 0.024483211 *
#> Dp_Dp 0.156306332 0.288355787 *
```

The inclusion probabilities of the $`\Pi`$ columns are one by
construction, as explained above, and the inclusion probabilities of the
remaining columns do not carry a period: variable selection decides on a
regressor for the whole sample.

### Coefficient and volatility paths

`plot` draws one figure per block of coefficients: the cointegration
matrix $`\Pi_t`$, the lagged differences, the deterministic terms and
the covariance matrix of the error term, each on its own.

``` r

plot(model)
```

![plot of chunk paths](figures/tvp-sv-vec-paths-1.png)

plot of chunk paths

![plot of chunk paths](figures/tvp-sv-vec-paths-2.png)

plot of chunk paths

![plot of chunk paths](figures/tvp-sv-vec-paths-3.png)

plot of chunk paths

![plot of chunk paths](figures/tvp-sv-vec-paths-4.png)

plot of chunk paths

In the figure of $`\Pi_t`$, reading down a column gives the response of
each differenced series to the level of one variable, and reading across
a row gives the error correction of one equation. A panel that sits at
zero among the $`\Gamma`$ coefficients is one that variable selection
has switched off.

A single element of $`\Pi_t`$ is easier to read on its own. The
following computes it from the draws of the loadings, which are the
first $`Kr`$ coefficients of each period of `posterior$a$coeffs`, and
the draws of the cointegration vectors in `posterior$beta$coeffs`.

``` r

pi_path <- function(object, equation, variable, ci = 0.9) {
  k <- object[["model"]][["k"]]
  r <- object[["model"]][["rank"]]
  m <- ncol(object[["data"]][["train"]][["z"]])
  k_ect <- ncol(object[["data"]][["train"]][["w"]])
  tt <- nrow(object[["data"]][["train"]][["y"]])
  draws <- nrow(object[["posterior"]][["a"]][["coeffs"]])

  i <- which(object[["model"]][["endogen"]] == equation)
  j <- which(dimnames(object[["data"]][["train"]][["w"]])[[2]] == variable)

  result <- matrix(NA_real_, draws, tt)
  for (period in 1:tt) {
    alpha <- object[["posterior"]][["a"]][["coeffs"]][, (period - 1) * m + 1:(k * r)]
    beta <- object[["posterior"]][["beta"]][["coeffs"]][, (period - 1) * k_ect * r + 1:(k_ect * r)]
    result[, period] <- rowSums(matrix(alpha, draws, k * r)[, k * 0:(r - 1) + i, drop = FALSE] *
                                  matrix(beta, draws, k_ect * r)[, k_ect * 0:(r - 1) + j, drop = FALSE])
  }

  bands <- t(apply(result, 2, stats::quantile,
                   probs = c((1 - ci) / 2, .5, 1 - (1 - ci) / 2)))

  stats::ts(bands, start = stats::start(object[["data"]][["train"]][["y"]]),
            frequency = stats::frequency(object[["data"]][["train"]][["y"]]))
}

path <- pi_path(model, equation = "r", variable = "l.r")

stats::ts.plot(path, lty = c(2, 1, 2),
               main = "Short-term rate: adjustment to its own level",
               ylab = "Element of Pi")
abline(h = 0, lty = "dotted")
```

![plot of chunk pi-path](figures/tvp-sv-vec-pi-path-1.png)

plot of chunk pi-path

The band of this element is widest in the early 1980s and narrows over
the sample, while its median rises from about -0.04 to zero: whatever
adjustment of the short-term rate to its own level there is fades out
over the four decades.

Note what is *not* plotted here. The path of $`\beta_t`$ on its own is
not a quantity to read: a cointegration vector is identified only up to
normalisation, and normalising it on a variable whose coefficient
wanders near zero produces a path that swings between large values
without anything having happened to the model. $`\Pi_t`$ carries the
same information and is identified.

The volatilities are extracted from the draws of the error precision
exactly as in the VAR case.

``` r

volatility <- function(object) {
  k <- object[["model"]][["k"]]
  kk <- k * k
  tt <- nrow(object[["data"]][["train"]][["y"]])
  draws <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]]

  result <- matrix(NA_real_, tt, k)
  for (i in 1:tt) {
    result[i, ] <- rowMeans(apply(draws[, (i - 1) * kk + 1:kk], 1,
                                  function(x) {sqrt(diag(solve(matrix(x, k))))}))
  }

  stats::ts(result, start = stats::start(object[["data"]][["train"]][["y"]]),
            frequency = stats::frequency(object[["data"]][["train"]][["y"]]),
            names = object[["model"]][["endogen"]])
}

plot(volatility(model),
     main = "Posterior mean of the residual standard deviations")
```

![plot of chunk volatility](figures/tvp-sv-vec-volatility-1.png)

plot of chunk volatility

For these data the volatility block does most of its work in the
equation of the short-term interest rate. Its residual standard
deviation falls from about a quarter of a percentage point in the early
1980s to a tenth of that once the rate approached its lower bound, with
a brief return during the financial crisis. A model that averaged over
the sample would report credible bands for the 2010s that are far too
wide. The volatilities of the other two equations move much less: the
one of inflation falls by about half until 2000, the one of the
long-term rate by about a quarter, most of it in the first half of the
1980s.

## A prior centred on the maximum likelihood estimate

The prior used above says nothing about where the cointegration space
lies: under $`\beta_t = \rho \beta_{t-1} + \eta_t`$ the marginal prior
of $`sp(\beta_t)`$ is uniform in every period. The working paper version
of Koop, León-González and Strachan (2011) shows how to centre it on a
given space $`sp(H)`$ instead, by putting
$`P_\tau = H H^{\prime} + H_\perp T H_\perp^{\prime}`$ into the
transition of the state equation,

``` math
\beta_t = \rho (I_r \otimes P_\tau) \beta_{t-1} + \eta_t, \qquad \eta_t \sim N(0, I).
```

The part of $`\beta_t`$ along $`sp(H)`$ keeps the autocorrelation
$`\rho`$ and the part off it decays faster, so the mode of the marginal
distribution of $`sp(\beta_t)`$ is $`sp(H)`$ in every period, while the
space of one period is still centred near the space of the period
before.

With `coint = list(rho = 0.999, p_tau_i = "ml", weight = 1)`,
`add_priors` takes $`H`$ from Johansen’s (1995) maximum likelihood
estimate of the space and chooses $`T`$ so that the prior spread of
$`sp(\beta_t)`$ around it in each period matches the sampling variance
of that estimate, with `weight` counting how many samples’ worth of
information the prior adds. The state before the sample gets the
stationary distribution of the new state equation. The estimate is
computed from the error correction term as it is stored in the model,
here on the scale of the data.

``` r

model_ml <- create_bvecmodel(data, p = 2, r = 1, const = "unrestricted",
                             tvp = TRUE, error = "sv+covar", varsel = "bvs",
                             iterations = 3000, burnin = 1000)

model_ml <- add_priors(model_ml,
                       coef = list(v_i = 1, v_i_det = 1 / 10,
                                   shape = 3, rate = 0.000001, rate_det = 0.01),
                       coint = list(rho = 0.999, p_tau_i = "ml", weight = 1),
                       sigma = list(shape = 3, rate = 0.01,
                                    mu = 0, v_i = 0.01,
                                    state_variance = 0.05, offset = 1e-8),
                       varsel = list(inprior = 0.5, exclude_det = TRUE))
```

How informative this prior can be is bounded by $`\rho`$: even $`T = 0`$
leaves the space $`1 - \rho^2`$ of room to move in each period. Where
the sampling variance of the estimate in some direction is smaller than
that, `add_priors` sets $`T`$ to zero in that direction and warns. On
these data both directions off $`sp(H)`$ stay above the bound, so the
prior matches the sampling variance without a warning. The eigenvalues
of the transition that was stored show the outcome – one for the
direction of $`sp(H)`$ itself, and one for each direction off it:

``` r

round(eigen(model_ml[["priors"]][["beta"]][["p_tau"]], symmetric = TRUE)$values, 3)
#> [1] 1.000 0.989 0.168
```

``` r

set.seed(20260913)
model_ml <- add_initial_values(model_ml)
model_ml <- add_posterior_coefficients(model_ml)
```

The same element of $`\Pi_t`$ as above, with the bands under the
noninformative prior in grey:

``` r

path_ml <- pi_path(model_ml, equation = "r", variable = "l.r")

stats::ts.plot(path, path_ml, lty = c(2, 1, 2, 2, 1, 2),
               col = rep(c("grey60", "black"), each = 3),
               main = "Short-term rate: adjustment to its own level",
               ylab = "Element of Pi")
abline(h = 0, lty = "dotted")
legend("bottomleft", legend = c("Noninformative", "Centred on ML"),
       col = c("grey60", "black"), lty = 1, bty = "n")
```

![plot of chunk ml-prior-path](figures/tvp-sv-vec-ml-prior-path-1.png)

plot of chunk ml-prior-path

Averaged over the sample, the width of the 90% bands of every element of
$`\Pi_t`$ under the prior centred on the maximum likelihood estimate,
relative to the noninformative one – rows are the variables in the error
correction term, columns the equations:

``` r

band_width <- function(object) {
  sapply(object[["model"]][["endogen"]], function(equation) {
    sapply(dimnames(object[["data"]][["train"]][["w"]])[[2]], function(variable) {
      bands <- pi_path(object, equation = equation, variable = variable)
      mean(bands[, 3] - bands[, 1])
    })
  })
}

round(band_width(model_ml) / band_width(model), 2)
#>         r   lr   Dp
#> l.r  0.76 0.72 0.86
#> l.lr 0.83 0.72 1.00
#> l.Dp 0.82 0.82 1.08
```

A ratio below one means that the prior centred on the maximum likelihood
estimate narrows the band of that element. Here it does so for all
elements but two, most in the equation of the long-term rate. The band
of the adjustment of inflation to the lagged long-term rate is
unchanged, and that of inflation to its own level is about eight percent
wider than under the noninformative prior. It need not narrow every
band, and that is a property of how the prior is built rather than a
failure of it: a larger `weight` makes the prior tighter by lowering
$`T`$, but $`T`$ is also how much of the position of the space carries
over from one period to the next. In a direction where $`T`$ is zero,
the position does not carry over at all, each period is informed by its
own data alone, and the posterior in that direction can be wider than
under the noninformative prior – see
[`?add_priors.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md).
How this plays out for the elements of $`\Pi_t`$ depends on how the
directions of the space line up with the variables.

Two further things to keep in mind. First, the centre is estimated from
the same sample the model is estimated on, so the data enter twice, and
where the bands do narrow they overstate how precisely the space is
known. Second, and most important for these data, $`H`$ is the estimate
of a *constant* space. Centring every period on it pulls the path
towards the average of the sample, which works against the question the
model is estimated to answer, whether the long-run relationship moved. A
prior like this is best suited to a model whose cointegration space is
expected to move around a stable centre, rather than to drift away from
one.

## Impulse responses and forecasts

Impulse responses and variance decompositions of a VEC model are
obtained from its VAR representation in levels, which `vec_to_var`
computes draw by draw. A time varying model is a VEC model per period,
so the transformation is applied period by period, and the result is a
VAR model with time varying parameters:

``` r

bvar_form <- vec_to_var(model)

bvar_form[["model"]][["algorithm"]]
#> [1] "VarTvpStochvol"
```

The draws of the variances of the state equations of the VEC
coefficients and of $`\rho`$ are not carried over: they describe how the
coefficients of the VEC model drift and have no counterpart among the
coefficients in levels. The paths of the coefficients and the draws of
the error term are, and they are all the application functions need.

Impulse responses and forecast error variance decompositions are those
of one period, which argument `period` selects, the last period being
the default. The response of inflation to a shock to the short-term
interest rate in 1985Q1:

``` r

feir <- irf(bvar_form, impulse = "r", response = "Dp", n_ahead = 20,
            period = period_of(model, 1985, 1))

plot(feir, main = "Forecast error impulse response, 1985Q1",
     xlab = "Period", ylab = "Response")
```

![plot of chunk feir-1985](figures/tvp-sv-vec-feir-1985-1.png)

plot of chunk feir-1985

Forecasts are made from the VEC model itself. `add_forecast_input`
prepares the regressors of the forecast periods, in levels, and
`add_posterior_forecasts` simulates each draw forward from the last
period of the sample: the loadings, the short-run coefficients and the
volatilities take a step of their random walks every quarter, the
cointegration vectors a step of their state equation, and the VAR in
levels is rebuilt from them in every period. The bands therefore reflect
both the uncertainty about where the coefficients are at the end of the
sample and the drift the model allows after it.

``` r

model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)

plot(predict(model, n_ahead = 8), n_pre = 20)
```

![plot of chunk forecasts](figures/tvp-sv-vec-forecasts-1.png)

plot of chunk forecasts

![plot of chunk forecasts](figures/tvp-sv-vec-forecasts-2.png)

plot of chunk forecasts

![plot of chunk forecasts](figures/tvp-sv-vec-forecasts-3.png)

plot of chunk forecasts

The VAR representation above has no state equation for its coefficients
in levels, so it can only hold them at their values in the last period.
`forecast_states = "hold"` gives that forecast from the VEC model as
well. The width of the 90% intervals eight quarters ahead shows what the
drift adds:

``` r

held <- add_posterior_forecasts(model, forecast_states = "hold")

interval_width <- function(object) {
  fcst <- predict(object, n_ahead = 8)[["fcst"]]
  apply(fcst[8, , ], 1, function(x) diff(stats::quantile(x, c(0.05, 0.95))))
}

round(rbind(simulated = interval_width(model), held = interval_width(held)), 2)
#>              r   lr   Dp
#> simulated 1.82 2.40 1.69
#> held      0.61 1.19 1.07
```

Eight quarters ahead, carrying the drift forward about doubles the
interval of the long-term rate, widens that of inflation by half and
makes that of the short-term rate between two and three times as wide. A
forecast that holds the coefficients at the end of the sample
understates the uncertainty this model implies.

The constant coefficient model with stochastic volatility, estimated
with `tvp = FALSE`, remains the model a TVP-SV specification has to
beat. The two can be compared out of sample by estimating both over an
expanding window with `use_expanding_window` and forecasting the models
of its windows directly: `add_forecast_input` and
`add_posterior_forecasts` take VEC models, forecast them in levels and,
for the time varying one, simulate its cointegration vectors and
volatility forward over the forecast horizon.

## Citing bvartools

If you use `bvartools` in published work, please cite it.
`citation("bvartools")` prints the reference, and the package has the
DOI [10.5281/zenodo.22736604](https://doi.org/10.5281/zenodo.22736604),
which always resolves to the latest archived version.

## References

Johansen, S. (1995). *Likelihood-based inference in cointegrated vector
autoregressive models*. Oxford: Oxford University Press.

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference
in a time varying cointegration model. *Journal of Econometrics,
165*(2), 210-220. <https://doi.org/10.1016/j.jeconom.2011.07.007>

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224-242.
<https://doi.org/10.1080/07474930903382208>

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204-230.
<https://doi.org/10.1002/jae.1271>

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Primiceri, G. E. (2005). Time varying structural vector autoregressions
and monetary policy. *Review of Economic Studies, 72*(3), 821-852.
<https://doi.org/10.1111/j.1467-937X.2005.00353.x>
