# Impulse Responses and Variance Decompositions of a Stored Model

Impulse response functions and forecast error variance decompositions of
a model in an HDF5 file, computed a piece of the chain at a time.

## Usage

``` r
# S3 method for class 'modelfile'
irf(x, ..., chunk = 100)

# S3 method for class 'modelfile'
fevd(x, ..., chunk = 100)
```

## Arguments

- x:

  an object of class 'modelfile', from
  [`read_model_file`](https://rdrr.io/pkg/bayests/man/read_model_file.html).

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

For a variance decomposition the pieces hand back the shares of each of
their draws, and the decomposition is summarised once over all of them,
so `statistic`, `ci`, `normalise_gir` and `max_groups` describe the
whole chain rather than the piece they were applied to.

An identification given per draw is read at the draw's position in the
whole chain: a list of impact matrices is cut to each piece's draws, and
an impact function is called with the draw's number in the chain. A
piece in which no draw is identified – under `type = "sign"`, or an
impact function that declines all of its draws – contributes nothing.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
[`add_predictive_loglik.expandingwindow()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.expandingwindow.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`pool_forecasts()`](https://franzmohr.github.io/bvartools/reference/pool_forecasts.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md)

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
stored <- read_model_file(file)

ir <- irf(stored, impulse = "invest", response = "income", n_ahead = 5,
          chunk = 25)
shares <- fevd(stored, response = "income", n_ahead = 5, chunk = 25)
```
