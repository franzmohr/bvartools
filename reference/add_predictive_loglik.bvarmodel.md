# Add the Log Predictive Density of a Forecast

Scores the forecast of an object of class 'bvarmodel' against the
observations its horizon realised.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_predictive_loglik(object, test_sample = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel', usually, the result of a call to
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md).

- test_sample:

  a time-series object used as test data. If `NULL` (default), the
  values in `data$test$y` of the object are used, which is what a model
  carries after it has been scored once and what a model read from a
  file was written with.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with `posterior$forecast$loglik` added, a
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) object with one row per
draw and one column per scored period, and with the values it was scored
against in `data$test$y`.

## Details

The draws of the log predictive density of period \\T + i\\ are stored
in `posterior$forecast$loglik`, one row per draw and one column per
scored period, beside the paths in `posterior$forecast$forecasts` that
they score. The numbers are BayesTS's own: this method fills
`data$test$y` and hands the object to the same C++ that computes the
score when the `bayests` programme is run on a model file.

Each column conditions on the observations the periods before it
realised, not on the path the forecast simulated. The log of the mean of
the draws of a column is therefore the one step ahead predictive density
given everything known up to that period, and those sum over the columns
to \\\ln p(y\_{T+1}, \ldots, y\_{T+h} \| y\_{1}, \ldots, y\_{T})\\, the
log predictive likelihood of the whole realised stretch, which
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
reports as `LPL`. Scoring against a simulated history instead would give
the marginal density of each period on its own, which is a different
quantity and one whose columns could not be added up.

With the history realised rather than simulated, the regressors of the
scored periods do not depend on the draw. The score is then the model's
own pointwise log-likelihood over those periods – the same expression as
[`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md)
evaluates, over a different sample – so nothing is written down twice
and the two cannot drift apart. What still moves with the draw is
everything that follows a state equation: time varying coefficients,
time varying error covariances and stochastic volatilities take one step
of their random walk per scored period, in the order and subject to the
same `model$forecast_states` that
[`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
steps them in. Under `"simulate"` those steps are drawn, so the score is
drawn too and two runs give two answers, exactly as two runs of a
forecast do.

Fewer realised periods than the horizon is not an error. The periods
that are there are the ones that can be scored, and the rest of the
forecast is left alone. A `test_sample` that does not reach the forecast
at all leaves the object unchanged, which is what a model estimated to
the end of a series looks like.

Structural models are refused: their regressors include the
contemporaneous observations, so the realised row is not built from its
lags alone, and their density carries the Jacobian of \\A_0\\ besides.
So are the asymmetric Laplace algorithms of quantile estimation, which
do not forecast in the first place.

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
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
data("e1")
orig <- diff(log(e1)) * 100
train <- window(orig, end = c(1978, 4))

# Create model
model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

model <- add_initial_values(model)
model <- add_posterior_coefficients(model)

# Forecast the four periods that were held back
model <- add_forecast_input(model, n_ahead = 4)
model <- add_posterior_forecasts(model)

# Score them against what those periods realised
model <- add_predictive_loglik(model, test_sample = orig)
dim(model[["posterior"]][["forecast"]][["loglik"]])
#> [1] 20  4

# The log predictive likelihood is criterion "LPL"
selection_criteria(model)[["LPL"]]
#>        mean median    qlower   qupper
#> 1 -25.05172     NA -31.19764 -18.9058
```
