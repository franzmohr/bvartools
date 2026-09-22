# Impulse Responses and Variance Decompositions of a Stored Model

Impulse response functions and forecast error variance decompositions of
a model in an HDF5 file, computed a piece of the chain at a time.

## Usage

``` r
# S3 method for class 'bvarfile'
irf(x, ..., chunk = 100)

# S3 method for class 'bvarfile'
fevd(x, ..., chunk = 100)
```

## Arguments

- x:

  an object of class 'bvarfile', from
  [`open_model`](https://franzmohr.github.io/bvartools/reference/open_model.md).

- ...:

  arguments of
  [`irf.bvarmodel`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md)
  or
  [`fevd.bvarmodel`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
  which do the work on each piece.

- chunk:

  how many draws are read at a time. Defaults to 100.

## Value

What
[`irf.bvarmodel`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md)
and
[`fevd.bvarmodel`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md)
return: an object of class 'bvarirf' or 'bvarfevd'.

## Details

Both quantities are sums over the draws: an impulse response is a set of
quantiles of the responses of the draws, a variance decomposition their
mean. Neither needs the draws together, so the chain is read in pieces
of `chunk`, each piece contributes what it has, and the pieces are put
together at the end. What comes out is what the same call on the model
in memory gives.

For a variance decomposition the pieces are combined before the shares
are normalised and before groups are collapsed, so `normalise_gir` and
`max_groups` describe the whole chain rather than the piece they were
applied to.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
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

data("e1")
e1 <- diff(log(e1)) * 100
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
model <- add_posterior_coefficients(add_initial_values(model))

file <- file.path(tempdir(), "bvartools-example-irf.h5")
unlink(file)
write_to_hdf5(model, filename = file)
stored <- open_model(file)

ir <- irf(stored, impulse = "invest", response = "income", n_ahead = 5,
          chunk = 25)
shares <- fevd(stored, response = "income", n_ahead = 5, chunk = 25)
```
