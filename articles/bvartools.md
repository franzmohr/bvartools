# Introduction to bvartools

## Introduction

The package `bvartools` implements functions for Bayesian inference of
linear time series models with a focus on vector autoregressive (VAR)
and vector error correction (VEC). It separates a typical analytics
workflow into multiple steps:

- *Model set-up*
  - Produce data matrices for given lag orders and model types, which
    can be used for posterior simulation. This includes a set-up for
    expanding window estimation.
  - Set prior hyperparameters either manually or based on established
    approaches such as the Minnesota prior.
  - Set initial values based on maximum likelihood estimates or drawing
    from a prior.
- *Posterior simulation*
  - Perform posterior simulation of coefficients, forecasts and
    likelihoods using algorithms from the
    [BayesTS](https://github.com/franzmohr/BayesTS) library.
  - Researchers can also choose to use they own algorithms.
- *Evaluation*
  - Traditional summary statistics for individual coefficients
  - In-sample performance measures (LL, AIC, BIC, HQ, WAIC, LOOIC)
  - Out-of-sample performance (MAFE, RMSFE)
- *Application*
  - Forecasts
  - Impulse response functions
  - Forecast error variance decomposition

In each step researchers are provided with the opportunity to fine-tune
a model according to their specific requirements or to use the default
framework for commonly used models and priors. Since version 1.0.0 the
package includes simulation functions from the
[BayesTS](https://github.com/franzmohr/BayesTS) C++ library.

For Bayesian inference of *VAR models* the package covers

- Standard BVAR models with independent normal-Wishart priors
- BVAR models employing stochastic search variable selection à la
  Gerorge, Sun and Ni (2008)
- BVAR models employing Bayesian variable selection à la Korobilis
  (2013)
- Structural BVAR models, where the structural coefficients are
  estimated from contemporary endogenous variables (A-model)
- Stochastic volatility (SV) of the errors à la Kim, Shephard and Chip
  (1998)
- Time varying parameter models (TVP-VAR)
- Bayesian quantile VARs, which estimate a conditional quantile instead
  of the conditional mean, à la Kozumi and Kobayashi (2011)

For Bayesian inference of *cointegrated VAR models* the package
implements the algorithm of Koop, León-González and Strachan (2010)
\[KLS\] – which places identification restrictions on the cointegration
space – in the following variants

- The BVEC model as presented in Koop, León-González and Strachan (2010)
- The KLS model employing stochastic search variable selection à la
  Gerorge, Sun and Ni (2008)
- The KLS modol employing Bayesian variable selection à la Korobilis
  (2013)
- Structural BVEC models, where the structural coefficients are
  estimated from contemporaneous endogenous variables (A-model).
  However, no further restrictions are made regarding the cointegration
  term.
- Stochastic volatility (SV) of the errors à la Kim, Shephard and Chip
  (1998)
- Time varying parameter models (TVP-VEC) à la Koop, León-González and
  Strachan (2011)\[^tvpvec\]

For Bayesian inference of *dynamic factor models* the package implements
the althorithm used in the textbook of Chan, Koop, Poirer and Tobias
(2019).

This introduction to `bvartools` provides the code to set up and
estimate a basic Bayesian VAR (BVAR) model.[^1] The first part covers a
basic workflow, where the standard posterior simulation algorithm of the
package is employed for Bayesian inference. The second part presents a
workflow for a posterior algorithm as it could be implemented by a
researcher.

## Model set-up

### Data

For both illustrations the data set E1 from Lütkepohl (2006) is used. It
contains data on West German fixed investment, disposable income and
consumption expenditures in billions of DM from 1960Q1 to 1982Q4. Like
in the textbook only the first 73 observations of the log-differenced
series are used.

``` r

library(bvartools)

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Reduce number of oberservations
e1 <- window(e1, end = c(1978, 4))

# Plot the series
plot(e1)
```

![plot of chunk data](figures/bvartools-data-1.png)

plot of chunk data

### Setting up a model

The `create_bvarmodel` function produces an object, which contains
information on the specification of the VAR model that should be
estimated. The following code specifies a VAR(2) model with an intercept
term. The number of iterations and burn-in draws is already specified at
this stage.

``` r

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 5000, burnin = 1000)
```

Note that the function is also capable of generating more than one
model. For example, specifying `p = 0:2` would result in three models.

### Adding prior hyperparameters

Function `add_priors` produces priors for the specified model(s) in
object `model` and augments the object accordingly.

``` r

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
```

If researchers want to fine-tune individual prior specifications, this
can be done by directly accessing the respective elements in object
`model_with_prior`.

### Adding initial values

Function `add_initial_values` adds initial values of posterior
coefficients to `model_with_prior`. The default behaviour of the
function is to obtain an LS estimate of the coefficients.

``` r

# The seed of the posterior simulation is drawn here
set.seed(1234567)
model <- add_initial_values(model, method = "ols")
```

## Obtaining posterior draws

`bvartools` comes with built-in posterior simulation functions from the
[BayesTS](https://github.com/franzmohr/BayesTS) library. By using
function `add_posterior_coefficients` the model input is forwarded to a
posterior function and the output is added to the original object:

``` r

model <- add_posterior_coefficients(model)
```

If researchers prefer to use their own posterior algorithms, this can be
done by specifying argument `posterior_function` with a function that
uses object `model` as its input.

## Evaluation

### Inspect posterior draws

Posterior draws can be visually inspected by using the `plot` function.
By default, it produces a series of histograms of all estimated
coefficients.

``` r

plot(model)
```

![plot of chunk histograms](figures/bvartools-histograms-1.png)

plot of chunk histograms

![plot of chunk histograms](figures/bvartools-histograms-2.png)

plot of chunk histograms

![plot of chunk histograms](figures/bvartools-histograms-3.png)

plot of chunk histograms

Alternatively, the trace plot of the post-burnin draws can be drawn by
adding the argument `type = "trace"`:

``` r

plot(model, type = "trace")
```

![plot of chunk traceplots](figures/bvartools-traceplots-1.png)

plot of chunk traceplots

![plot of chunk traceplots](figures/bvartools-traceplots-2.png)

plot of chunk traceplots

![plot of chunk traceplots](figures/bvartools-traceplots-3.png)

plot of chunk traceplots

### Summary statistics

Summary statistics can be obtained in the usual way using the `summary`
method.

``` r

summary(model)
#> 
#> Bayesian VAR model with p = 2 
#> 
#> Endogenous variables: invest, income, cons
#> 
#> Variable: invest 
#> 
#>              Mean     SD Naive SD Time-series SD    2.5%     50%    97.5%  
#> invest.l1 -0.3203 0.1280 0.001810       0.001810 -0.5719 -0.3201 -0.07193 *
#> income.l1  0.1424 0.5594 0.007912       0.007108 -0.9692  0.1442  1.23388  
#> cons.l1    0.9642 0.6725 0.009511       0.009646 -0.3367  0.9506  2.31280  
#> invest.l2 -0.1614 0.1295 0.001832       0.001776 -0.4122 -0.1613  0.09020  
#> income.l2  0.1086 0.5532 0.007823       0.007823 -0.9933  0.1097  1.21134  
#> cons.l2    0.9361 0.6707 0.009486       0.009739 -0.3763  0.9410  2.25209  
#> const     -1.6550 1.7451 0.024679       0.024679 -5.1629 -1.6257  1.67947  
#> 
#> Variable: income 
#> 
#>                Mean      SD  Naive SD Time-series SD     2.5%      50%  97.5%  
#> invest.l1  0.043943 0.03298 0.0004665      0.0004665 -0.01952  0.04372 0.1098  
#> income.l1 -0.152671 0.14073 0.0019903      0.0019903 -0.43313 -0.15387 0.1204  
#> cons.l1    0.287758 0.17147 0.0024249      0.0023351 -0.04064  0.28651 0.6269  
#> invest.l2  0.049701 0.03254 0.0004602      0.0004602 -0.01377  0.04985 0.1144  
#> income.l2  0.019545 0.13897 0.0019653      0.0019653 -0.25122  0.01942 0.2960  
#> cons.l2   -0.009856 0.17233 0.0024371      0.0024371 -0.34588 -0.01156 0.3257  
#> const      1.579222 0.45619 0.0064515      0.0064515  0.69652  1.57580 2.4841 *
#> 
#> Variable: cons 
#> 
#>                Mean      SD  Naive SD Time-series SD      2.5%       50%
#> invest.l1 -0.001837 0.02677 0.0003786      0.0003882 -0.054492 -0.001704
#> income.l1  0.228460 0.11544 0.0016325      0.0014185  0.001675  0.228818
#> cons.l1   -0.268581 0.13948 0.0019726      0.0019726 -0.543926 -0.266272
#> invest.l2  0.033545 0.02635 0.0003726      0.0003726 -0.017927  0.033485
#> income.l2  0.356565 0.11183 0.0015814      0.0015814  0.138655  0.356735
#> cons.l2   -0.025501 0.14022 0.0019831      0.0019831 -0.295589 -0.027811
#> const      1.300183 0.36488 0.0051602      0.0051602  0.585607  1.296486
#>               97.5%  
#> invest.l1  0.049108  
#> income.l1  0.450706 *
#> cons.l1   -0.001353 *
#> invest.l2  0.085652  
#> income.l2  0.573624 *
#> cons.l2    0.250518  
#> const      2.028044 *
#> 
#> Variance-covariance matrix:
#> 
#>                  Mean     SD Naive SD Time-series SD    2.5%     50%  97.5%  
#> invest_invest 22.3351 4.0218 0.056876       0.062536 15.7824 21.9124 31.176 *
#> invest_income  0.7227 0.7380 0.010438       0.011472 -0.7273  0.7060  2.235  
#> invest_cons    1.2908 0.6069 0.008583       0.009741  0.1795  1.2550  2.615 *
#> income_income  1.4446 0.2589 0.003662       0.004044  1.0184  1.4159  2.032 *
#> income_cons    0.6434 0.1685 0.002383       0.002631  0.3615  0.6277  1.009 *
#> cons_cons      0.9353 0.1683 0.002381       0.002624  0.6632  0.9153  1.316 *
```

As expected for an algrotihm with uninformative priors the posterior
means are fairly close to the results of the frequentist estimator,
which can be obtaind in the following way:

``` r

# Obtain data for LS estimator
k <- model[["model"]][["k"]]
y <- t(model[["data"]][["train"]][["y"]])
x <- t(model[["data"]][["train"]][["x"]])

# Calculate LS estimates
A_freq <- tcrossprod(y, x) %*% solve(tcrossprod(x))

# Round estimates and print
round(A_freq, 3)
#>        invest.l1 income.l1 cons.l1 invest.l2 income.l2 cons.l2  const
#> invest    -0.320     0.146   0.961    -0.161     0.115   0.934 -1.672
#> income     0.044    -0.153   0.289     0.050     0.019  -0.010  1.577
#> cons      -0.002     0.225  -0.264     0.034     0.355  -0.022  1.293
```

### Thin results

The MCMC series in object `model` can be thinned using

``` r

model <- thin(model, thin = 10)
```

## Application

### Forecasts

Forecasts are obtained in two steps. First, function
`add_forecast_input` generates the data of the forecast periods, where
deterministic terms can be provided in the argument `deterministic` and
unmodelled variables in the argument `exogen`. If they are not provided,
the function tries to obtain them from the model and gives a message, if
it does not succeed. Second, function `add_posterior_forecasts` then
simulates the forecasts.

``` r

model <- add_forecast_input(model, n_ahead = 10)
model <- add_posterior_forecasts(model)
```

Forecasts with credible bands can be extracted from the model using the
`predict` function. The `n_ahead` argument can be used to set the given
forecast horizon. However, the number of available forecast periods is
limited to specification used in `add_forecast_input`.

``` r

bvar_pred <- predict(model, n_ahead = 10)

plot(bvar_pred, n_pre = 20)
```

![plot of chunk forecasts](figures/bvartools-forecasts-1.png)

plot of chunk forecasts

![plot of chunk forecasts](figures/bvartools-forecasts-2.png)

plot of chunk forecasts

![plot of chunk forecasts](figures/bvartools-forecasts-3.png)

plot of chunk forecasts

### Impulse response analysis

`bvartools` supports commonly used impulse response functions. See
<https://www.r-econometrics.com/timeseries/irf/> for an introduction.

#### Forecast error impulse response

``` r

FEIR <- irf(model, impulse = "income", response = "cons", n_ahead = 8)

plot(FEIR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
```

![plot of chunk feir](figures/bvartools-feir-1.png)

plot of chunk feir

#### Orthogonalised impulse response

``` r

OIR <- irf(model, impulse = "income", response = "cons", n_ahead = 8, type = "oir")

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response")
```

![plot of chunk oir](figures/bvartools-oir-1.png)

plot of chunk oir

#### Generalised impulse response

``` r

GIR <- irf(model, impulse = "income", response = "cons", n_ahead = 8, type = "gir")

plot(GIR, main = "Generalised Impulse Response", xlab = "Period", ylab = "Response")
```

![plot of chunk gir](figures/bvartools-gir-1.png)

plot of chunk gir

### Variance decomposition

`bvartools` also supports forecast error variance decomposition (FEVD)
and generalised forecast error variance decomposition.

#### Forecast error variance decomposition

``` r

bvar_fevd_oir <- fevd(model, response = "cons")

plot(bvar_fevd_oir, main = "OIR-based FEVD of consumption")
```

![plot of chunk fevd-oir](figures/bvartools-fevd-oir-1.png)

plot of chunk fevd-oir

#### Generalised forecast error variance decomposition

It is also possible to calculate FEVDs, which are based on generalised
impulse responses (GIR). Note that these do not automatically add up to
unity. However, this could be changed by adding `normalise_gir = TRUE`
to the function’s arguments.

``` r

bvar_fevd_gir <- fevd(model, response = "cons", type = "gir")

plot(bvar_fevd_gir, main = "GIR-based FEVD of consumption")
```

![plot of chunk fevd-gir](figures/bvartools-fevd-gir-1.png)

plot of chunk fevd-gir

## Using `bvartools` with user-written algorithms

`bvartools` was created to assist researchers in building and evaluating
their own posterior simulation algorithms for linear BVAR models.
Functions `create_bvarmodel` simply help to quickly obtain the relevant
data matrices for posterior simulation. Estimation can be done using
algortihms, which are usually implemented by researchers themselves. But
once posterior draws are obtained `bvartools` can assist throughout in
the following steps of the analysis. In this context the main
contributions of the package are:

- Function `bvar` collects the output of a Gibbs sampler in a
  standardized object, which can be used for subsequent steps in an
  analysis.
- Functions such as `predict`, `irf`, `fevd` for forecasting, impulse
  response analysis and forecast error variance decomposition,
  respectively, use the output of `bvar` and, hence, researchers do not
  have to implement them themselves and can save time.
- Computationally intensive functions - such as for posterior
  simulation - are written in C++ using the `RcppArmadillo` package of
  Eddelbuettel and Sanderson (2014).[^2] This decreases calculation time
  and makes the code less complex and, thus, less prone to mistakes.

If researchers are willing to rely on the model generation and
evaluation functions of `bvartools`, the only remaining step is to
illustrate how user specific algorithms can be combined with the
functional framework of the package. This is shown in the remainder of
this introduction.

### Model set-up and prior specifications

These steps are exactly the same as described above. Thus, the following
Gibbs sampler departs from object `model_with_priors` from above.

### A Gibbs sampler algorithm

``` r

# Reset random number generator for reproducibility
set.seed(1234567)

iterations <- 10000 # Number of saved iterations of the Gibbs sampler
burnin <- 5000 # Number of burn-in draws
draws <- iterations + burnin # Total number of MCMC draws

y <- t(model[["data"]][["train"]][["y"]])
x <- t(model[["data"]][["train"]][["x"]])

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
m <- k * nrow(x) # Number of estimated coefficients

# Set (uninformative) priors
a_mu_prior <- model[["priors"]][["a"]][["mu"]] # Vector of prior parameter means
a_v_i_prior <- model[["priors"]][["a"]][["v_inv"]] # Inverse of the prior covariance matrix

u_sigma_df_prior <- model[["priors"]][["u_sigma"]][["df"]] # Prior degrees of freedom
u_sigma_scale_prior <- model[["priors"]][["u_sigma"]][["scale"]] # Prior covariance matrix
u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom

# Initial values
u_sigma_i <- diag(1 / .00001, k)

# Data containers for posterior draws
draws_a <- matrix(NA, m, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  # Draw conditional mean parameters
  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
  
  # Draw variance-covariance matrix
  u <- y - matrix(a, k) %*% x # Obtain residuals
  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
  
  # Store draws
  if (draw > burnin) {
    draws_a[, draw - burnin] <- a
    # Store inverted Sigma_i to obtain Sigma
    draws_sigma[, draw - burnin] <- solve(u_sigma_i)
  }
}
```

### `bvarmodel` objects

The `bvar` function can be used to collect relevant output of the Gibbs
sampler into a standardized object of class `bvarmodel`, which can be
used by functions such as `predict` to obtain forecasts or `irf` for
impulse respons analysis.

``` r

y <- model[["data"]][["train"]][["y"]]

bvar_est_two <- bvar(y = y,
                     x = model[["data"]][["train"]][["x"]],
                     A = draws_a[1:18,],
                     C = draws_a[19:21, ],
                     Sigma = draws_sigma)
```

Since the output of
[`add_posterior_coefficients()`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
is an object of class `bvarmodel` as well, the calculation of summary
statistics, forecasts, impulse responses and forecast error variance
decompositions is performed as described above.

``` r

summary(bvar_est_two)
#> 
#> Bayesian VAR model with p = 2 
#> 
#> Endogenous variables: invest, income, cons
#> 
#> Variable: invest 
#> 
#>              Mean     SD Naive SD Time-series SD    2.5%     50%    97.5%  
#> invest.l1 -0.3210 0.1280 0.001280       0.001286 -0.5734 -0.3216 -0.06831 *
#> income.l1  0.1472 0.5654 0.005654       0.005654 -0.9602  0.1439  1.26083  
#> cons.l1    0.9662 0.6768 0.006768       0.006778 -0.3453  0.9560  2.32034  
#> invest.l2 -0.1600 0.1263 0.001263       0.001305 -0.4049 -0.1612  0.08747  
#> income.l2  0.1036 0.5519 0.005519       0.005439 -0.9653  0.1010  1.19676  
#> cons.l2    0.9348 0.6894 0.006894       0.006894 -0.4165  0.9283  2.28239  
#> const     -1.6637 1.7556 0.017556       0.017556 -5.1131 -1.6428  1.81572  
#> 
#> Variable: income 
#> 
#>                Mean      SD  Naive SD Time-series SD     2.5%       50%  97.5%
#> invest.l1  0.043539 0.03283 0.0003283      0.0003340 -0.02134  0.043640 0.1078
#> income.l1 -0.152587 0.14272 0.0014272      0.0014354 -0.43484 -0.152102 0.1278
#> cons.l1    0.287003 0.17215 0.0017215      0.0017583 -0.05264  0.284572 0.6301
#> invest.l2  0.049836 0.03215 0.0003215      0.0003215 -0.01315  0.049738 0.1135
#> income.l2  0.019209 0.13846 0.0013846      0.0013846 -0.25074  0.020273 0.2888
#> cons.l2   -0.008994 0.17079 0.0017079      0.0017079 -0.34237 -0.009633 0.3335
#> const      1.577324 0.44978 0.0044978      0.0044978  0.69837  1.573286 2.4624
#>            
#> invest.l1  
#> income.l1  
#> cons.l1    
#> invest.l2  
#> income.l2  
#> cons.l2    
#> const     *
#> 
#> Variable: cons 
#> 
#>                Mean      SD  Naive SD Time-series SD      2.5%       50%
#> invest.l1 -0.002623 0.02648 0.0002648      0.0002648 -0.054699 -0.002433
#> income.l1  0.223178 0.11668 0.0011668      0.0011668 -0.003841  0.222297
#> cons.l1   -0.263006 0.13888 0.0013888      0.0013888 -0.539179 -0.262530
#> invest.l2  0.033789 0.02612 0.0002612      0.0002612 -0.017709  0.033990
#> income.l2  0.354398 0.11138 0.0011138      0.0011302  0.131559  0.356159
#> cons.l2   -0.020351 0.13878 0.0013878      0.0013661 -0.294508 -0.019763
#> const      1.292296 0.35786 0.0035786      0.0035786  0.590719  1.290012
#>              97.5%  
#> invest.l1 0.049704  
#> income.l1 0.449267  
#> cons.l1   0.007515  
#> invest.l2 0.085769  
#> income.l2 0.571058 *
#> cons.l2   0.254939  
#> const     2.006569 *
#> 
#> Variance-covariance matrix:
#> 
#>                  Mean     SD Naive SD Time-series SD    2.5%     50%  97.5%  
#> invest_invest 22.3072 4.0178 0.040178       0.044990 15.8399 21.7999 31.327 *
#> invest_income  0.7561 0.7313 0.007313       0.008046 -0.6330  0.7286  2.294  
#> invest_cons    1.2956 0.6022 0.006022       0.006781  0.2157  1.2701  2.576 *
#> income_income  1.4378 0.2610 0.002610       0.002854  1.0191  1.4057  2.038 *
#> income_cons    0.6442 0.1690 0.001690       0.001873  0.3596  0.6280  1.032 *
#> cons_cons      0.9348 0.1681 0.001681       0.001872  0.6610  0.9141  1.312 *
```

## Citing bvartools

If you use `bvartools` in published work, please cite it.
`citation("bvartools")` prints the reference, and the package has the
DOI [10.5281/zenodo.22736604](https://doi.org/10.5281/zenodo.22736604),
which always resolves to the latest archived version.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). *Bayesian
Econometric Methods* (2nd ed.). Cambridge: University Press.

Eddelbuettel, D., & Sanderson C. (2014). RcppArmadillo: Accelerating R
with high-performance C++ linear algebra. *Computational Statistics and
Data Analysis, 71*, 1054-1063.
<https://doi.org/10.1016/j.csda.2013.02.005>

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553-580.
<https://doi.org/10.1016/j.jeconom.2007.08.017>

Kim, S., Shephard, N., & Chib, S. (1998). Stochastic volatility:
Likelihood inference and comparison with ARCH models. *Review of
Economic Studies, 65*(3), 361-393.
<https://doi.org/10.1111/1467-937X.00050>

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224-242.
<https://doi.org/10.1080/07474930903382208>

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference
in a time varying cointegration model. *Journal of Econometrics,
165*(2), 210-220. <https://doi.org/10.1016/j.jeconom.2011.07.007>

Koop, G., Pesaran, M. H., & Potter, S.M. (1996). Impulse response
analysis in nonlinear multivariate models. *Journal of Econometrics
74*(1), 119-147. <https://doi.org/10.1016/0304-4076(95)01753-4>

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204-230.
<https://doi.org/10.1002/jae.1271>

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis
in linear multivariate models. *Economics Letters, 58*(1), 17-29.
<https://doi.org/10.1016/S0165-1765(97)00214-0>

Sanderson, C., & Curtin, R. (2016). Armadillo: a template-based C++
library for linear algebra. *Journal of Open Source Software, 1*(2), 26.
<https://doi.org/10.21105/joss.00026>

[^1]: Further examples about the use of the `bvartools` package are
    available at <https://www.r-econometrics.com/timeseriesintro/>.

[^2]: `RcppArmadillo` is the `Rcpp` bridge to the open source
    ‘Armadillo’ library of Sanderson and Curtin (2016).
