# Open a Folder of Stored Models

Opens a folder written with
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
and returns a handle to the models in it, which the estimation steps
work on one at a time rather than reading all of them into the session.

## Usage

``` r
open_models(folder)

# S3 method for class 'bvarfolder'
print(x, ...)

model_files(x, models = NULL)
```

## Arguments

- folder:

  path to a folder holding models in HDF5 files.

- x:

  an object of class 'bvarfolder'.

- ...:

  further arguments passed to or from other methods.

- models:

  character vector of the models to name, as the manifest names them, or
  `NULL` for all of them.

## Value

A list of class 'bvarfolder' with the elements `folder` and `manifest`,
a data frame with the columns `model`, `file`, `group` and the
specification of each model.

## Details

[`read_models_from_folder`](https://franzmohr.github.io/bvartools/reference/read_models_from_folder.md)
reads a folder into a list, draws and all, which is what to do with
models that fit. Many do not. A lag and rank grid is dozens of models,
an expanding window is one per quarter, and a model with time varying
coefficients keeps a coefficient path per draw, so a folder can be tens
of gigabytes while each model in it is a few hundred megabytes.

A handle carries what the folder holds – a row per model, saying where
it is and what it is – and none of the draws.
[`map_models`](https://franzmohr.github.io/bvartools/reference/map_models.md)
and the estimation steps for a handle then read one model, change it,
write it back and drop it, so a step costs one model per worker rather
than the whole folder.
[`open_model`](https://franzmohr.github.io/bvartools/reference/open_model.md)
is the same idea for one model whose chain is too long to hold: a handle
to a file, read a piece at a time.

The manifest is built by reading the specification of every model, which
is a set of attributes rather than any of its data, so opening a folder
is cheap. A model is named by where its file sits below `folder`,
without the extension, and by its group where a file holds more than one
– the names
[`read_models_from_folder`](https://franzmohr.github.io/bvartools/reference/read_models_from_folder.md)
gives.

`model_files` names the files the models lie in, which is what a program
that works on the files is pointed at. With
[`bayests_files`](https://franzmohr.github.io/bvartools/reference/bayests_files.md)
a folder is estimated as

    run <- bayests_files(executable = "/opt/bayests/bin/bayests")
    run(model_files(stored), jobs = 6)

and the draws never pass through R. Files are named once however many
models they hold, since BayesTS works through every model in a file it
is given.

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
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

data("e1")
train <- diff(log(e1)) * 100

models <- create_bvarmodel(data = train, p = 1:3, deterministic = "const",
                           iterations = 20, burnin = 10)

folder <- file.path(tempdir(), "bvartools-example-folder")
unlink(folder, recursive = TRUE)
dir.create(folder, recursive = TRUE)
write_to_hdf5(models, folder = folder)

stored <- open_models(folder)
stored
#> 3 models in /tmp/RtmpcU1dnC/bvartools-example-folder 
#>   3 VarNormalWishart 
stored[["manifest"]][, c("model", "p", "iterations")]
#>                                       model p iterations
#> 1 VarNormalWishart-varsel=none-p=01-n=1-001 1         20
#> 2 VarNormalWishart-varsel=none-p=02-n=1-001 2         20
#> 3 VarNormalWishart-varsel=none-p=03-n=1-001 3         20
```
