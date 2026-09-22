# Add Predictive Log-Likelihood

Adds the draws of the one-step-ahead log predictive density to the
windows of an expanding window exercise, which
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
sums to the log predictive likelihood.

## Usage

``` r
add_predictive_loglik(object, ...)

# S3 method for class 'expandingwindow'
add_predictive_loglik(object, ...)

# S3 method for class 'modellist'
add_predictive_loglik(object, ...)
```

## Arguments

- object:

  an object of class `"expandingwindow"` whose windows hold posterior
  draws, or a `"modellist"` of such objects, as returned by
  [`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md).

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object`, with element `predictive` added to every window
but the last. It is a list with `loglik`, the draws of the log
predictive density, and `period`, the time of the predicted observation.

## Details

The windows of an expanding window exercise grow by one period at a
time, so window \\i + 1\\ holds exactly one observation that window
\\i\\ has not seen. For every window but the last, the function
evaluates the density of that observation, \\\Delta y\_{t}\\, given the
regressors of period \\t\\, at each posterior draw of the window, and
stores the draws in element `predictive` of the window, together with
the period \\t\\. The log of their mean is the log predictive density
\\\ln p(\Delta y_t \| y\_{t-1}, \ldots, y_1)\\, and its sum over the
windows is the log predictive likelihood of Geweke and Amisano (2011),
which Koop, León-González and Strachan (2011) use to choose between time
varying cointegration models and their ranks.

A draw of a model with constant coefficients and a constant error
covariance describes period \\t\\ as it describes the sample. Everything
that follows a state equation is carried one period forward with that
equation first:

- time varying coefficients by one step of their random walk with the
  drawn state variance, and a coefficient that variable selection
  excluded stays at zero;

- a time varying cointegration space by \\\beta_t = \rho (I_r \otimes
  P\_\tau) \beta\_{t-1} + \eta_t\\, \\\eta_t \sim N(0, I)\\, with the
  drawn \\\rho\\ or the one of the prior, where \\P\_\tau\\ is the
  identity unless the prior centres the space;

- time varying error covariances by one step of their random walk;

- stochastic volatilities by one step of the random walk of the log
  variances. Its state variance is taken from
  `posterior$u_sigma_inv$sigma` if the sampler stored it, and is
  otherwise drawn from its conditional posterior given the drawn path of
  the log variances and the prior in `priors$u_sigma`.

Unlike the pointwise log-likelihood of
[`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md),
whose states have seen the observation they are evaluated at, the
predictive density conditions only on the data before period \\t\\. This
is what makes it the criterion for models whose coefficients or
variances follow a state equation, where leave-one-out importance
sampling fails for exactly the periods a path bends towards.

The function is available for VEC models, which includes the rank zero
models in differences that ranks are compared with, and for VAR models,
whose regressors are the lags of the endogenous variables, the exogenous
variables and the deterministic terms of the period that is predicted. A
VAR model is the case of a VEC model without an error correction term,
so the density is the same expression with the cointegration block left
out. The windows of a VEC model must have been simulated on error
correction terms that are neither scaled nor centred, or put back with
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
first, and structural models are not supported.

The expression is the normal density of the observation, so the function
is available for the algorithms whose observation is normal given the
parameters of a draw, and refuses the others. The asymmetric Laplace
algorithms of quantile estimation, `"VarNormalAld"` and `"VarTvpAld"`,
are refused: the precision their samplers store is the one of the normal
that their scale mixture conditions on period by period, not the density
of an observation, whose mixing variable would have to be integrated
out.

The methods for a single fitted model,
[`add_predictive_loglik.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md)
and
[`add_predictive_loglik.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
score a forecast against the observations its horizon realised instead.
That is the same statistic over a different set of periods – one step
ahead densities, each conditioning on the realised history before it –
and it is BayesTS that computes it rather than the R code here. It is
stored in `posterior$forecast$loglik`, beside the forecasts it scores,
rather than in `predictive`.

## References

Geweke, J., & Amisano, G. (2011). Hierarchical Markov normal mixture
models with applications to financial asset returns. *Journal of Applied
Econometrics, 26*(1), 1–29.
[doi:10.1002/jae.1119](https://doi.org/10.1002/jae.1119)

Koop, G., León-González, R., & Strachan, R. W. (2011). Bayesian
inference in a time varying cointegration model. *Journal of
Econometrics, 165*(2), 210–220.
[doi:10.1016/j.jeconom.2011.07.007](https://doi.org/10.1016/j.jeconom.2011.07.007)

## See also

Methods for a single fitted model:
[`add_predictive_loglik.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md).

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`analysis_of_stored_models`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`folder_steps`](https://franzmohr.github.io/bvartools/reference/folder_steps.md),
[`map_draws()`](https://franzmohr.github.io/bvartools/reference/map_draws.md),
[`map_models()`](https://franzmohr.github.io/bvartools/reference/map_models.md),
[`open_model()`](https://franzmohr.github.io/bvartools/reference/open_model.md),
[`open_models()`](https://franzmohr.github.io/bvartools/reference/open_models.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data("e6")
e6 <- e6 * 100

# Create model
model <- create_bvecmodel(e6, p = 2, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 0.0001))

# Estimate the model on the last four expanding windows
model <- use_expanding_window(model, start = c(1998, 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

# One-step-ahead log predictive densities
model <- add_predictive_loglik(model)

# The log predictive likelihood is criterion "LPL"
selection_criteria(model)
#> 
#> 
#> ------------------------------------------
#> Predictive
#> ------------------------------------------
#> 
#>  Criterion   Mean Quantile (2.5%) Quantile (97.5%)
#>        LPL -9.138          -10.84           -7.433
#> 
#> One-step-ahead log predictive likelihood over 4 periods; numerical standard error 0.08726.
```
