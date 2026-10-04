# Pool Forecasts

Combines the forecasts of several models into one equal-weight pool,
which is compared with the models like any of them.

## Usage

``` r
pool_forecasts(...)

# S3 method for class 'forecastpool'
add_priors(object, ...)

# S3 method for class 'forecastpool'
add_initial_values(object, ...)

# S3 method for class 'forecastpool'
add_seed(object, seed, ...)

# S3 method for class 'forecastpool'
add_posterior_coefficients(object, ...)

# S3 method for class 'forecastpool'
add_posterior_loglik(object, ...)

# S3 method for class 'forecastpool'
add_forecast_input(object, ...)

# S3 method for class 'forecastpool'
add_posterior_forecasts(object, ...)

# S3 method for class 'forecastpool'
add_predictive_loglik(object, ...)

# S3 method for class 'forecastpool'
thin(x, thin = 10, ...)

# S3 method for class 'poolwindow'
add_priors(object, ...)

# S3 method for class 'poolwindow'
add_initial_values(object, ...)

# S3 method for class 'poolwindow'
add_seed(object, seed, ...)

# S3 method for class 'poolwindow'
add_posterior_coefficients(object, ...)

# S3 method for class 'poolwindow'
add_posterior_loglik(object, ...)

# S3 method for class 'poolwindow'
add_forecast_input(object, ...)

# S3 method for class 'poolwindow'
add_posterior_forecasts(object, ...)

# S3 method for class 'poolwindow'
add_predictive_loglik(object, ...)

# S3 method for class 'poolwindow'
thin(x, thin = 10, ...)

# S3 method for class 'forecastpool'
write_to_hdf5(object, ...)

# S3 method for class 'poolwindow'
write_to_hdf5(object, ...)

# S3 method for class 'forecastpool'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'poolwindow'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- ...:

  two or more models with forecasts: objects of class
  `"expandingwindow"`, or single models of class `"bvarmodel"` or
  `"bvecmodel"`, or `"modellist"` objects of either. Arguments may be
  named; the names label the members of the pool.

- object, x:

  an object of class `"forecastpool"` or `"poolwindow"`.

- seed:

  not used, since a pool is not simulated.

- thin:

  an integer specifying the thinning interval between successive draws.

- digits:

  the number of significant digits.

## Value

An object of class `"forecastpool"`, which behaves like an expanding
window exercise, or, if single models were pooled, a single pooled
forecast of class `"poolwindow"`. Each pooled forecast holds the draws
of the pool in `posterior$forecast$forecasts` and, in `model`, the names
of the members in `members` and the number of draws taken from each in
`draws`.

## Details

A pool of forecasts is a mixture of predictive distributions: the
predictive distribution of the pool is the average of those of its
members. Combining forecasts in this way is one of the most reliable
ways of improving them, because the errors of different models are not
perfectly correlated and no single specification is best in every period
(Bates and Granger, 1969; Hall and Mitchell, 2007; Geweke and Amisano,
2011).

The pool is built from draws the members already hold. For every
forecast, the same number of draws is taken at random from each member –
as many as the member with the fewest draws has – and stacked, so that
every member has the same weight. The draws are selected with R's random
number generator, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
makes a pool reproducible.

**The members must forecast the same thing.** The pool contains the
endogenous variables all members share, and the forecast horizons all of
them reach. A variable forecast by only some of the members is left out,
and a pool of models that share no variable is refused. A VEC model
forecasts the levels of its variables, so it can be pooled with VAR
models of the same levels, but not with VAR models of their growth
rates. Forecasts aggregated to annual figures can only be pooled with
forecasts aggregated in the same way.

The windows of expanding window exercises are matched by the end of
their estimation samples, so that each pooled forecast combines
forecasts made with the same information. Windows that not every member
has are left out, with a message. Single models and expanding windows
cannot be pooled with each other, and external forecasts, which are
point forecasts rather than predictive distributions, cannot be pooled
at all: they are compared beside a pool, by combining both with
[`combine_models`](https://franzmohr.github.io/bvartools/reference/combine_models.md).

The log predictive densities that
[`add_predictive_loglik`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md)
and the forecasts of BayesTS provide are pooled as the forecasts are,
which gives the log of the average predictive density of the members –
the score of the mixture, not of its draws. They are kept only if every
member carries them and all members have the same endogenous variables,
because a predictive density is a joint density of all the variables of
a model.

The weights are equal. Weights that favour members which forecast better
in the past must only use the periods before each forecast, or they
flatter the pool with information it did not have.

**A pool is formed after the members' forecasts are drawn**, and the
functions that estimate a model leave it unchanged. It can be passed to
[`add_forecast_errors`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
and
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md),
and included in a list of models with
[`combine_models`](https://franzmohr.github.io/bvartools/reference/combine_models.md).
It cannot be written to a file: write its members and form the pool
again after reading them back.

A pool is formed from draws its members already hold, so the functions
that add priors, initial values, posterior draws, forecasts or
predictive log-likelihoods to a model return it unchanged. This allows
to combine a pool with its members in a list of class `"modellist"` and
to apply the usual workflow to all elements of that list.

## References

Bates, J. M., & Granger, C. W. J. (1969). The combination of forecasts.
*Journal of the Operational Research Society, 20*(4), 451–468.
[doi:10.1057/jors.1969.103](https://doi.org/10.1057/jors.1969.103)

Geweke, J., & Amisano, G. (2011). Optimal prediction pools. *Journal of
Econometrics, 164*(1), 130–141.
[doi:10.1016/j.jeconom.2011.02.017](https://doi.org/10.1016/j.jeconom.2011.02.017)

Hall, S. G., & Mitchell, J. (2007). Combining density forecasts.
*International Journal of Forecasting, 23*(1), 1–13.
[doi:10.1016/j.ijforecast.2006.08.001](https://doi.org/10.1016/j.ijforecast.2006.08.001)

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
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
data("e1")
e1 <- diff(log(e1)) * 100
train <- window(e1, end = c(1978, 4))

# Two specifications, each estimated over expanding windows
model <- create_bvarmodel(train, p = 1:2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 1 / 10, v_i_det = 1 / 100),
                    sigma = list(df = "k", scale = 1))
model <- use_expanding_window(model, start = c(1978, 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model)
model <- add_forecast_input(model, n_ahead = 2)
model <- add_posterior_forecasts(model)

# The equal-weight pool of both
set.seed(123)
pool <- pool_forecasts(p1 = model[[1]], p2 = model[[2]])

# Scored and compared like the models
all_models <- combine_models(model, pool)
all_models <- add_forecast_errors(all_models, test_sample = e1)
selection_criteria(all_models)
#> 
#> 
#> ------------------------------------------
#> Out-of-sample
#> ------------------------------------------
#> 
#> Mean absolute forecast errors (MAFE)
#> 
#>  Variable h Model 1 Model 2 Model 3
#>    invest 1   3.662  4.0380   3.850
#>    income 1   1.143  1.1740   1.159
#>      cons 1   1.302  0.9985   1.150
#>    invest 2   4.363  4.9460   4.655
#>    income 2   1.343  1.1702   1.257
#>      cons 2   1.146  1.1920   1.169
#> 
#> 
#> Root mean squared forecast errors (RMSFE)
#> 
#>  Variable h Model 1 Model 2 Model 3
#>    invest 1   4.914   4.992   4.953
#>    income 1   1.394   1.427   1.411
#>      cons 1   1.560   1.254   1.415
#>    invest 2   5.558   6.563   6.081
#>    income 2   1.613   1.495   1.555
#>      cons 2   1.440   1.466   1.454
```
