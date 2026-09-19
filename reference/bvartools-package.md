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
'bvarmodel' or 'bvecmodel' with the elements

- `data`:

  the data matrices, with the estimation sample in `data$train` (`y`,
  `x` and `z` in SUR form).

- `model`:

  the specification: variables, lag orders, deterministic terms, error
  and variable selection types, iterations and burn-in draws.

- `priors`:

  the prior hyperparameters, added by
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

- `initial`:

  the starting values of the sampler, added by
  [`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md).

- `posterior`:

  the draws added by
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  and the later `add_posterior_*` functions.

Posterior draws are stored as
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) objects with one row
per draw and one column per parameter.
[`bvar`](https://franzmohr.github.io/bvartools/reference/bvar.md) and
[`bvec`](https://franzmohr.github.io/bvartools/reference/bvec.md), which
collect the draws of a sampler written by the user, expect the
transpose: one row per parameter and one column per draw.

Passing a vector to an argument such as `p` or `r` creates one model per
specification in a list of class 'modellist', and
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
creates a list of class 'expandingwindow'. The functions of the workflow
accept these lists as well and apply to each model in turn.

## Workflow of a VAR model

1.  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
    builds the data matrices from a time-series object. Lag order,
    deterministic terms, exogenous variables and the number of
    iterations and burn-in draws are set here.
    [`transform_variables`](https://franzmohr.github.io/bvartools/reference/transform_variables.md)
    applies the transformation codes of FRED-MD and FRED-QD beforehand.

2.  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
    adds prior hyperparameters. Its arguments `coef` and `sigma` are
    named lists without defaults, and neither may be empty: `coef` needs
    `v_i` or `minnesota`, and `sigma` the elements of the chosen
    `error`, which are `df` and `scale` for `"wishart"` and `shape` and
    `rate` for `"gamma"` and `"ald"`. Structural models cannot use the
    Wishart prior.

3.  [`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
    adds starting values, by default from a least squares estimate.

4.  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
    runs the sampler.
    [`plot`](https://franzmohr.github.io/bvartools/reference/plot.bvarmodel.md),
    [`summary`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md)
    and
    [`thin`](https://franzmohr.github.io/bvartools/reference/thin.bvarmodel.md)
    inspect and thin the draws.

5.  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md)
    followed by
    [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
    simulates forecasts, which
    [`predict`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
    summarises.
    [`irf`](https://franzmohr.github.io/bvartools/reference/irf.md),
    [`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md)
    and
    [`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md)
    compute impulse responses, variance decompositions and connectedness
    measures.

## Choosing a model

The type of model is fixed in
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
and
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
and the priors that fit it in
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md):

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
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

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
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).
Forecasts, impulse responses, variance decompositions, spillovers and
sign restrictions are computed for a 'bvarmodel', so an estimated VEC
model is first converted to its VAR representation in levels with
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md).

## Identification

[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md) and
[`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md)
provide forecast error, orthogonalised and generalised impulse responses
through their `type` argument.
[`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md)
identifies the shocks of an estimated model by the signs of their
impulse responses instead.

## Model comparison and forecast evaluation

- In sample:
  [`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md),
  then
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  (LL, AIC, BIC, HQ, WAIC, LOOIC), then
  [`choose_best_model`](https://franzmohr.github.io/bvartools/reference/choose_best_model.md).

- Out of sample:
  [`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
  before the priors are added, the usual estimation and forecasting
  steps, then
  [`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
  with the test sample and
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  for MAFE and RMSFE.
  [`plot_forecast_errors_by_period`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md)
  plots the errors over time.

- [`combine_models`](https://franzmohr.github.io/bvartools/reference/combine_models.md)
  and
  [`align_model_obs`](https://franzmohr.github.io/bvartools/reference/align_model_obs.md)
  put differently specified models on a common sample, and
  [`create_external_forecast`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  brings forecasts produced elsewhere into the same comparison.

## User-written samplers

[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
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

[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
and
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
save and restore a model with its draws, one per file or several side by
side in groups listed by
[`list_models_in_hdf5`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md).
The files follow the format of the BayesTS command line, so a model can
also be estimated there and read back.
[`read_models_from_folder`](https://franzmohr.github.io/bvartools/reference/read_models_from_folder.md)
and
[`read_expanding_window_model_from_folder`](https://franzmohr.github.io/bvartools/reference/read_expanding_window_model_from_folder.md)
read a whole folder.

## Common mistakes

- Calling a step without assigning its result, which discards what it
  added.

- Calling
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
  without
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md).
  The `n_ahead` given there is also the longest horizon `predict`
  returns; a longer one is shortened to it.

- Calling
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  for in-sample criteria without
  [`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md).

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

- Report bugs at <https://github.com/franzmohr/bvartools/issues>

## Author

**Maintainer**: Franz X. Mohr <franz.x.mohr@outlook.com>
([ORCID](https://orcid.org/0009-0003-8890-7781))

Authors:

- Franz X. Mohr <franz.x.mohr@outlook.com>
  ([ORCID](https://orcid.org/0009-0003-8890-7781))

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
#>                  Mean        SD   Naive SD Time-series SD       2.5%
#> invest.l1 -0.26289162 0.1016472 0.01437508     0.01437508 -0.4594323
#> income.l1  0.24843939 0.5992033 0.08474014     0.06587082 -0.9530594
#> cons.l1    0.59589536 0.6240210 0.08824990     0.05035741 -0.6442137
#> invest.l2 -0.09746837 0.1151917 0.01629056     0.02224465 -0.3034178
#> income.l2  0.13181345 0.4917103 0.06953833     0.06953833 -0.6698627
#> cons.l2    0.52929174 0.5805670 0.08210458     0.08210458 -0.6357997
#> const     -0.54169520 1.5238423 0.21550384     0.21550384 -3.8419171
#>                   50%      97.5%  
#> invest.l1 -0.26671431 -0.0945254 *
#> income.l1  0.23569900  1.4312909  
#> cons.l1    0.66884428  1.4920327  
#> invest.l2 -0.11287320  0.1094395  
#> income.l2  0.09122699  1.0271221  
#> cons.l2    0.62485237  1.4374361  
#> const     -0.51809318  2.1081975  
#> 
#> Variable: income 
#> 
#>                  Mean         SD    Naive SD Time-series SD         2.5%
#> invest.l1  0.04426586 0.02773175 0.003921862    0.002896493 -0.005609147
#> income.l1 -0.13698966 0.12547322 0.017744594    0.017744594 -0.363279191
#> cons.l1    0.29358676 0.13147666 0.018593607    0.018593607  0.041933253
#> invest.l2  0.05624536 0.03419211 0.004835494    0.005980819 -0.006199024
#> income.l2  0.02187446 0.12972204 0.018345466    0.018345466 -0.206694065
#> cons.l2    0.02982776 0.14695360 0.020782378    0.020782378 -0.215527248
#> const      1.39958846 0.32042655 0.045315158    0.045315158  0.820276297
#>                   50%      97.5%  
#> invest.l1  0.05025868 0.09093887  
#> income.l1 -0.15232899 0.12998444  
#> cons.l1    0.30924916 0.51615217 *
#> invest.l2  0.05327211 0.11596999  
#> income.l2  0.01338268 0.21731427  
#> cons.l2    0.02386611 0.31361871  
#> const      1.40618306 1.97001463 *
#> 
#> Variable: cons 
#> 
#>                   Mean         SD    Naive SD Time-series SD         2.5%
#> invest.l1  0.006208935 0.02699946 0.003818301    0.003818301 -0.054861083
#> income.l1  0.291576211 0.10710157 0.015146449    0.015146449  0.100980356
#> cons.l1   -0.289857993 0.12100249 0.017112337    0.017112337 -0.505685347
#> invest.l2  0.043452315 0.02418062 0.003419656    0.004369468 -0.002115442
#> income.l2  0.367576715 0.10469687 0.014806373    0.014806373  0.188687914
#> cons.l2   -0.118673536 0.13169306 0.018624210    0.018624210 -0.348117325
#> const      1.296523183 0.30982827 0.043816334    0.043816334  0.706984389
#>                    50%       97.5%  
#> invest.l1  0.008528816  0.04718163  
#> income.l1  0.299026172  0.52428883 *
#> cons.l1   -0.303202288 -0.06223235 *
#> invest.l2  0.042450683  0.08636827  
#> income.l2  0.358325761  0.54257735 *
#> cons.l2   -0.119002560  0.11526141  
#> const      1.287265993  1.86558869 *
#> 
#> Variance-covariance matrix:
#> 
#>                     Mean        SD   Naive SD Time-series SD       2.5%
#> invest_invest 20.4340370 3.4673940 0.49036356     0.60575744 15.3622913
#> invest_income  0.6754984 0.5833800 0.08250239     0.08250239 -0.2503046
#> invest_cons    1.3464135 0.5285908 0.07475402     0.07475402  0.5142435
#> income_income  1.3404008 0.2280484 0.03225092     0.03225092  0.9781828
#> income_cons    0.6648747 0.1706583 0.02413472     0.03400581  0.4105721
#> cons_cons      0.9883308 0.1800381 0.02546124     0.02546124  0.7437618
#>                      50%     97.5%  
#> invest_invest 19.8889678 28.194741 *
#> invest_income  0.5719998  1.982575  
#> invest_cons    1.2697250  2.623203 *
#> income_income  1.3111596  1.810006 *
#> income_cons    0.6520988  1.054431 *
#> cons_cons      0.9979511  1.406629 *
#> 

# Forecasts and impulse responses
model <- add_forecast_input(model, n_ahead = 8)
model <- add_posterior_forecasts(model)
pred <- predict(model, n_ahead = 8)
oir <- irf(model, impulse = "income", response = "cons", n_ahead = 8, type = "oir")
```
