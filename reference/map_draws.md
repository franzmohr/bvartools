# Apply a Function to the Draws of a Stored Model

Reads the chain of a stored model in pieces and applies a function to
the model holding each piece, which is how the analysis of a model too
large to hold is done.

## Usage

``` r
map_draws(x, f, ..., chunk = 100)
```

## Arguments

- x:

  an object of class 'bvarfile', from
  [`open_model`](https://franzmohr.github.io/bvartools/reference/open_model.md).

- f:

  a function taking a model and returning whatever the caller wants to
  keep of that piece of the chain.

- ...:

  further arguments passed to `f`.

- chunk:

  how many draws are read at a time. Defaults to 100.

## Value

A list with one element per piece of the chain.

## Details

The model `f` is given is the model in the file with the draws of one
piece, so anything that works on a model works on it. What it returns is
collected in a list, one element per piece, and combining those is the
caller's business: the responses of an impulse response are stacked, the
shares of a variance decomposition are averaged over the pieces they
came from.

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
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
model <- add_posterior_coefficients(add_initial_values(model))

file <- file.path(tempdir(), "bvartools-example-draws.h5")
unlink(file)
write_to_hdf5(model, filename = file)

stored <- open_model(file)

# The mean of every coefficient, without reading the chain at once
sums <- map_draws(stored, function(model) {
  colSums(unclass(model[["posterior"]][["a"]][["coeffs"]]))
}, chunk = 25)
means <- Reduce(`+`, sums) / stored[["draws"]]
```
