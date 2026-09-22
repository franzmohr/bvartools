# Estimation Steps on a Folder of Models

The steps of an estimation, applied to the models of a folder one at a
time rather than to a list of models in the session. Each reads a model,
changes it, writes it back and drops it, so a step costs one model per
worker.

## Usage

``` r
# S3 method for class 'bvarfolder'
add_priors(object, ..., cores = 1)

# S3 method for class 'bvarfolder'
add_initial_values(object, ..., cores = 1)

# S3 method for class 'bvarfolder'
add_seed(object, seed, ..., cores = 1)

# S3 method for class 'bvarfolder'
add_posterior_coefficients(object, ..., cores = 1, export = NULL)

# S3 method for class 'bvarfolder'
add_posterior_loglik(object, ..., cores = 1)

# S3 method for class 'bvarfolder'
thin(x, ..., cores = 1)

# S3 method for class 'bvarfolder'
selection_criteria(object, ..., cores = 1)
```

## Arguments

- object:

  an object of class 'bvarfolder', from
  [`open_models`](https://franzmohr.github.io/bvartools/reference/open_models.md).

- ...:

  further arguments passed to the method for a single model.

- cores:

  the number of worker processes. Defaults to one.

- seed:

  an integer, the seed of the first model. The models are numbered
  through in the order of the manifest, as the models of a list are.

- export:

  character vector of names of objects the workers need, for instance
  the sampler that `posterior_function` refers to. See
  [`map_models`](https://franzmohr.github.io/bvartools/reference/map_models.md).

- x:

  an object of class 'bvarfolder'.

## Value

The object in `object`, with its manifest brought up to date.

## Details

The methods do what their counterparts for a 'bvarmodel' or a
'bvecmodel' do, model by model, and they take the same arguments.

Each model is drawn with a seed of its own, so the draws do not depend
on how many workers a step runs on, nor on whether the models were
estimated in one run or in several. The seeds are the ones a list of
models gets from
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md):
`seed`, `seed + 1`, and so on in the order of the manifest.

A step takes `models`, which
[`map_models`](https://franzmohr.github.io/bvartools/reference/map_models.md)
understands: a character vector of the models it is applied to, the rest
of the folder being left as it is. An estimation of many expensive
models is run in parts that way, and one that was interrupted is picked
up where it stopped.

`add_posterior_coefficients` is the step to think twice about. It writes
the largest thing a model has, and a cluster of R workers carries every
draw into the session and out again to do it. Where the sampler is the
BayesTS executable,
[`bayests_files`](https://franzmohr.github.io/bvartools/reference/bayests_files.md)
runs it on the files instead:

    run <- bayests_files(executable = "/opt/bayests/bin/bayests")
    run(model_files(stored), jobs = 6)

which gives the same draws, since a model is drawn with the seed in its
file.

`selection_criteria` is the one method that reads rather than writes: it
hands back a criterion per model and leaves the folder as it is. It is
what
[`choose_best_model`](https://franzmohr.github.io/bvartools/reference/choose_best_model.md)
compares, so a grid of specifications too large to hold is still chosen
from.

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
train <- diff(log(e1)) * 100

models <- create_bvarmodel(data = train, p = 1:2, deterministic = "const",
                           iterations = 20, burnin = 10)

folder <- file.path(tempdir(), "bvartools-example-steps")
unlink(folder, recursive = TRUE)
dir.create(folder, recursive = TRUE)
write_to_hdf5(models, folder = folder)

stored <- open_models(folder)
stored <- add_priors(stored, coef = list(v_i = 1),
                     sigma = list(df = 3, scale = 1))
stored <- add_initial_values(stored)
stored <- add_seed(stored, 20260919)
stored <- add_posterior_coefficients(stored)
stored <- add_posterior_loglik(stored)
```
