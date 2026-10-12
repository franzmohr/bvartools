# bvartools: Bayesian Inference of Vector Autoregressive and Error Correction Models

Assists in the set-up of algorithms for Bayesian inference of vector
autoregressive (VAR) and error correction (VEC) models. Functions for
posterior simulation, forecasting, impulse response analysis and
forecast error variance decomposition are largely based on the
introductory texts of Chan, Koop, Poirier and Tobias (2019, ISBN:
9781108437493), Koop and Korobilis (2010)
[doi:10.1561/0800000013](https://doi.org/10.1561/0800000013) and
Luetkepohl (2006, ISBN: 9783540262398).

## Model objects

Every step of an analysis takes a model object and returns it with
something added, so the result of each call must be assigned back:
`model <- add_initial_values(model)`. A model object is a list of class
'bvarmodel' or 'bvecmodel' whose elements are `model`, the
specification, `data`, `priors`, `initial` and `posterior`, the draws –
stored as [`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) objects with
one row per draw and one column per parameter.

Passing a vector to an argument such as `p` or `r` creates one model per
specification in a list of class 'modellist', and
[`use_expanding_window`](https://rdrr.io/pkg/bayests/man/use_expanding_window.html)
creates a list of class 'expandingwindow'. The functions of the workflow
accept these lists as well and apply to each model in turn.

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the whole object element by element: what each one holds,
which step adds it, how the draws are laid out and which path of a model
file each corresponds to. It is the page to read before reaching into
the object by hand.

## Workflow of a VAR model

1.  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
    builds the data matrices from a time-series object. Lag order,
    deterministic terms, exogenous variables and the number of
    iterations and burn-in draws are set here.
    [`transform_variables`](https://franzmohr.github.io/bvartools/reference/transform_variables.md)
    applies the transformation codes of FRED-MD and FRED-QD beforehand.

2.  [`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html) adds
    prior hyperparameters. Its arguments `coef` and `sigma` are named
    lists without defaults, and neither may be empty: `coef` needs `v_i`
    or `minnesota`, and `sigma` the elements of the chosen `error`,
    which are `df` and `scale` for `"wishart"` and `shape` and `rate`
    for `"gamma"` and `"ald"`. Structural models cannot use the Wishart
    prior.

3.  [`add_initial_values`](https://rdrr.io/pkg/bayests/man/add_initial_values.html)
    adds starting values, by default from a least squares estimate.

4.  [`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
    runs the sampler.
    [`plot`](https://franzmohr.github.io/bvartools/reference/plot.bvarmodel.md),
    [`summary`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md)
    and
    [`thin`](https://franzmohr.github.io/bvartools/reference/thin.bvarmodel.md)
    inspect and thin the draws.

5.  [`add_forecast_input`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html)
    followed by
    [`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html)
    simulates forecasts, which
    [`predict`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
    summarises. [`irf`](https://rdrr.io/pkg/bayests/man/irf.html),
    [`fevd`](https://rdrr.io/pkg/bayests/man/fevd.html) and
    [`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md)
    compute impulse responses, variance decompositions and connectedness
    measures.

## Choosing a model

The type of model is fixed in
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
and
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
and the priors that fit it in
[`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html):

- `error` sets the error covariance: `"wishart"` (default), `"gamma"`
  for a diagonal covariance, `"sv"` for stochastic volatility and, for
  VAR models, `"ald"` for a quantile VAR whose quantiles are given in
  `quantile`. Quantile VARs estimate no error covariances and produce no
  forecasts.

- `tvp = TRUE` makes the coefficients time varying.

- `structural = TRUE` estimates the contemporaneous coefficients
  (A-model).

- `varsel` selects variables by `"ssvs"` (George et al., 2008) or
  `"bvs"` (Korobilis, 2013).

- A Minnesota prior is requested with `coef = list(minnesota = ...)` in
  [`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html).

[`minnesota_prior`](https://franzmohr.github.io/bvartools/reference/minnesota_prior.md),
[`ssvs_prior`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.md)
and
[`inclusion_prior`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.md)
return the same prior components on their own, for use in a sampler
written by the user.

## Error correction models

[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
additionally takes the cointegration rank `r` and whether constant,
trend and seasonal terms are `"restricted"` to the cointegration space
or `"unrestricted"`. Priors on the cointegration space go into the
`coint` argument of
[`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html).
Forecasts, impulse responses, variance decompositions, spillovers and
sign restrictions are computed for a 'bvarmodel', so an estimated VEC
model is first converted to its VAR representation in levels with
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md).

## Identification

[`irf`](https://rdrr.io/pkg/bayests/man/irf.html) and
[`fevd`](https://rdrr.io/pkg/bayests/man/fevd.html) provide forecast
error, orthogonalised and generalised impulse responses through their
`type` argument.
[`add_sign_restrictions`](https://rdrr.io/pkg/bayests/man/add_sign_restrictions.html)
identifies the shocks of an estimated model by the signs of their
impulse responses instead.

## Model comparison and forecast evaluation

- In sample:
  [`add_posterior_loglik`](https://rdrr.io/pkg/bayests/man/add_posterior_loglik.html),
  then
  [`selection_criteria`](https://rdrr.io/pkg/bayests/man/selection_criteria.html)
  (LL, AIC, BIC, HQ, WAIC, LOOIC), then
  [`choose_best_model`](https://rdrr.io/pkg/bayests/man/choose_best_model.html).

- Out of sample:
  [`use_expanding_window`](https://rdrr.io/pkg/bayests/man/use_expanding_window.html)
  before the priors are added, the usual estimation and forecasting
  steps, then
  [`add_forecast_errors`](https://rdrr.io/pkg/bayests/man/add_forecast_errors.html)
  with the test sample and
  [`selection_criteria`](https://rdrr.io/pkg/bayests/man/selection_criteria.html)
  for MAFE and RMSFE.
  [`plot_forecast_errors_by_period`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md)
  plots the errors over time.

- [`combine_models`](https://rdrr.io/pkg/bayests/man/combine_models.html)
  and
  [`align_model_obs`](https://rdrr.io/pkg/bayests/man/align_model_obs.html)
  put differently specified models on a common sample, and
  [`create_external_forecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  brings forecasts produced elsewhere into the same comparison.

## User-written samplers

[`add_posterior_coefficients`](https://rdrr.io/pkg/bayests/man/add_posterior_coefficients.html)
accepts a `posterior_function` that replaces the built-in sampler. For a
sampler written from scratch the package exports its building blocks,
among them
[`post_normal`](https://franzmohr.github.io/bvartools/reference/post_normal.md),
[`post_normal_sur`](https://franzmohr.github.io/bvartools/reference/post_normal_sur.md),
[`post_bvs`](https://franzmohr.github.io/bvartools/reference/post_bvs.md),
[`ssvs`](https://franzmohr.github.io/bvartools/reference/ssvs.md),
[`post_coint_kls`](https://franzmohr.github.io/bvartools/reference/post_coint_kls.md),
[`stochvol_ksc_1998`](https://franzmohr.github.io/bvartools/reference/stochvol_ksc_1998.md)
and
[`kalman_durbin_koopman_2002`](https://franzmohr.github.io/bvartools/reference/kalman_durbin_koopman_2002.md).
[`bvar`](https://franzmohr.github.io/bvartools/reference/bvar.md) and
[`bvec`](https://franzmohr.github.io/bvartools/reference/bvec.md)
collect the resulting draws in a model object, to which the rest of the
workflow applies.

## Storage

[`write_to_hdf5`](https://rdrr.io/pkg/bayests/man/write_to_hdf5.html)
and
[`read_model_from_hdf5`](https://rdrr.io/pkg/bayests/man/read_model_from_hdf5.html)
save and restore a model with its draws, one per file or several side by
side in groups listed by
[`read_model_groups`](https://rdrr.io/pkg/bayests/man/read_model_groups.html).
The files follow the format of the BayesTS command line, so a model can
also be estimated there and read back.
[`read_models_from_folder`](https://rdrr.io/pkg/bayests/man/read_models_from_folder.html)
and
[`read_expanding_window_model_from_folder`](https://rdrr.io/pkg/bayests/man/read_expanding_window_model_from_folder.html)
read a whole folder.

## Common mistakes

- Calling a step without assigning its result, which discards what it
  added.

- Calling
  [`add_posterior_forecasts`](https://rdrr.io/pkg/bayests/man/add_posterior_forecasts.html)
  without
  [`add_forecast_input`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html).
  The `n_ahead` given there is also the longest horizon `predict`
  returns; a longer one is shortened to it.

- Calling
  [`selection_criteria`](https://rdrr.io/pkg/bayests/man/selection_criteria.html)
  for in-sample criteria without
  [`add_posterior_loglik`](https://rdrr.io/pkg/bayests/man/add_posterior_loglik.html).

- Calling `predict`, `irf` or `fevd` on a 'bvecmodel' instead of the
  result of
  [`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md).

- Passing draws to
  [`bvar`](https://franzmohr.github.io/bvartools/reference/bvar.md) or
  [`bvec`](https://franzmohr.github.io/bvartools/reference/bvec.md) with
  one row per draw instead of one column per draw.

## Vignettes and data

Worked examples are in the vignettes, listed by
`browseVignettes("bvartools")`: `"bvartools"` (introduction),
`"minnesota-prior"`, `"ssvs"`, `"tvp-sv-var"`, `"quantile-var"`,
`"sign-restrictions"`, `"bvec"`, `"tvp-sv-vec"`, `"model-comparison"`
and `"horse-races"`. They use the data sets
[`e1`](https://franzmohr.github.io/bvartools/reference/e1.md),
[`e6`](https://franzmohr.github.io/bvartools/reference/e6.md),
[`us_macrodata`](https://franzmohr.github.io/bvartools/reference/us_macrodata.md)
and
[`uk_macrodata`](https://franzmohr.github.io/bvartools/reference/uk_macrodata.md).

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). *Bayesian
Econometric Methods* (2nd ed.). Cambridge: University Press.

Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation
smoother for state space time series analysis. *Biometrika, 89*(3),
603–615.

Eddelbuettel, D., & Sanderson C. (2014). RcppArmadillo: Accelerating R
with high-performance C++ linear algebra. *Computational Statistics and
Data Analysis, 71*, 1054–1063.
[doi:10.1016/j.csda.2013.02.005](https://doi.org/10.1016/j.csda.2013.02.005)

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Koop, G, & Korobilis, D. (2010), Bayesian multivariate time series
Methods for empirical macroeconomics, *Foundations and Trends in
Econometrics, 3*(4), 267–358.
[doi:10.1561/0800000013](https://doi.org/10.1561/0800000013)

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224–242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference
in a time varying cointegration model. *Journal of Econometrics,
165*(2), 210–220.
[doi:10.1016/j.jeconom.2011.07.007](https://doi.org/10.1016/j.jeconom.2011.07.007)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Primiceri, G. E. (2005). Time varying structural vector autoregressions
and monetary policy. *The Review of Economic Studies, 72*(3), 821–852.
[doi:10.1111/j.1467-937X.2005.00353.x](https://doi.org/10.1111/j.1467-937X.2005.00353.x)

Sanderson, C., & Curtin, R. (2016). Armadillo: a template-based C++
library for linear algebra. *Journal of Open Source Software, 1*(2), 26.
[doi:10.21105/joss.00026](https://doi.org/10.21105/joss.00026)

## See also

Useful links:

- <https://github.com/franzmohr/bvartools>

- <https://franzmohr.github.io/bvartools/>

- Report bugs at <https://github.com/franzmohr/bvartools/issues>

## Author

**Maintainer**: Franz X. Mohr <franz.x.mohr@outlook.com>
([ORCID](https://orcid.org/0009-0003-8890-7781)) \[copyright holder\]

Authors:

- Franz X. Mohr <franz.x.mohr@outlook.com>
  ([ORCID](https://orcid.org/0009-0003-8890-7781)) \[copyright holder\]

## Examples

``` r
data("e1")
e1 <- diff(log(e1)) * 100

# Set up, estimate and inspect a VAR(2) model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 50, burnin = 10)
# Number of iterations and burnin should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
summary(model)
#> 
#> Bayesian VAR model with p = 2 
#> 
#> Endogenous variables: invest, income, cons
#> 
#> Variable: invest 
#> 
#>               Mean     SD Naive SD Time-series SD    2.5%      50%    97.5%  
#> invest.l1 -0.26289 0.1016  0.01438        0.01438 -0.4594 -0.26671 -0.09453 *
#> income.l1  0.24844 0.5992  0.08474        0.06587 -0.9531  0.23570  1.43129  
#> cons.l1    0.59590 0.6240  0.08825        0.05036 -0.6442  0.66884  1.49203  
#> invest.l2 -0.09747 0.1152  0.01629        0.02224 -0.3034 -0.11287  0.10944  
#> income.l2  0.13181 0.4917  0.06954        0.06954 -0.6699  0.09123  1.02712  
#> cons.l2    0.52929 0.5806  0.08210        0.08210 -0.6358  0.62485  1.43744  
#> const     -0.54170 1.5238  0.21550        0.21550 -3.8419 -0.51809  2.10820  
#> 
#> Variable: income 
#> 
#>               Mean      SD Naive SD Time-series SD      2.5%      50%   97.5%  
#> invest.l1  0.04427 0.02773 0.003922       0.002896 -0.005609  0.05026 0.09094  
#> income.l1 -0.13699 0.12547 0.017745       0.017745 -0.363279 -0.15233 0.12998  
#> cons.l1    0.29359 0.13148 0.018594       0.018594  0.041933  0.30925 0.51615 *
#> invest.l2  0.05625 0.03419 0.004835       0.005981 -0.006199  0.05327 0.11597  
#> income.l2  0.02187 0.12972 0.018345       0.018345 -0.206694  0.01338 0.21731  
#> cons.l2    0.02983 0.14695 0.020782       0.020782 -0.215527  0.02387 0.31362  
#> const      1.39959 0.32043 0.045315       0.045315  0.820276  1.40618 1.97001 *
#> 
#> Variable: cons 
#> 
#>                Mean      SD Naive SD Time-series SD      2.5%       50%
#> invest.l1  0.006209 0.02700 0.003818       0.003818 -0.054861  0.008529
#> income.l1  0.291576 0.10710 0.015146       0.015146  0.100980  0.299026
#> cons.l1   -0.289858 0.12100 0.017112       0.017112 -0.505685 -0.303202
#> invest.l2  0.043452 0.02418 0.003420       0.004369 -0.002115  0.042451
#> income.l2  0.367577 0.10470 0.014806       0.014806  0.188688  0.358326
#> cons.l2   -0.118674 0.13169 0.018624       0.018624 -0.348117 -0.119003
#> const      1.296523 0.30983 0.043816       0.043816  0.706984  1.287266
#>              97.5%  
#> invest.l1  0.04718  
#> income.l1  0.52429 *
#> cons.l1   -0.06223 *
#> invest.l2  0.08637  
#> income.l2  0.54258 *
#> cons.l2    0.11526  
#> const      1.86559 *
#> 
#> Variance-covariance matrix:
#> 
#>                  Mean     SD Naive SD Time-series SD    2.5%     50%  97.5%  
#> invest_invest 20.4340 3.4674  0.49036        0.60576 15.3623 19.8890 28.195 *
#> invest_income  0.6755 0.5834  0.08250        0.08250 -0.2503  0.5720  1.983  
#> invest_cons    1.3464 0.5286  0.07475        0.07475  0.5142  1.2697  2.623 *
#> income_income  1.3404 0.2280  0.03225        0.03225  0.9782  1.3112  1.810 *
#> income_cons    0.6649 0.1707  0.02413        0.03401  0.4106  0.6521  1.054 *
#> cons_cons      0.9883 0.1800  0.02546        0.02546  0.7438  0.9980  1.407 *
#> 

# Forecasts and impulse responses
model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)
pred <- predict(model, n_ahead = 8)
oir <- irf(model, impulse = "income", response = "cons", n_ahead = 8, type = "oir")
```
