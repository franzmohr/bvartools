# Apply a Function to the Models of a Folder

Reads every model of a folder, applies a function to it, writes it back
and drops it, so that a step over many models costs the memory of one
model per worker rather than of all of them.

## Usage

``` r
map_models(x, f, ..., models = NULL, write = TRUE, cores = 1, export = NULL)
```

## Arguments

- x:

  an object of class 'bvarfolder', from
  [`open_models`](https://franzmohr.github.io/bvartools/reference/open_models.md).

- f:

  a function taking a model – a 'bvarmodel' or a 'bvecmodel' – and
  returning one. If it has an argument named `index`, the position of
  the model in the manifest is passed to it, which is what numbers seeds
  and what lets a caller tell the models apart.

- ...:

  further arguments passed to `f`.

- models:

  character vector of the models the step is applied to, as the manifest
  names them, or `NULL`, the default, for all of them. The others are
  left untouched, and so are their rows of the manifest.

- write:

  whether what `f` returns is written back. `TRUE`, the default, is a
  step of an estimation: `f` takes a model and returns one, and the
  model in the file becomes what it returned. `FALSE` reads only: `f`
  may return anything, nothing is written, and the results are collected
  and returned.

- cores:

  the number of worker processes the models are handled on. Defaults to
  one. With more than one the models are spread over a socket cluster;
  every model has its own file, so the workers share nothing.

- export:

  character vector of names of objects in the calling environment that
  the workers need, for instance the sampler that `posterior_function`
  refers to. Ignored with one core.

## Value

With `write = TRUE`, the object in `x` with its manifest brought up to
date, invisibly. With `write = FALSE`, a named list of what `f`
returned, one element per model.

## Details

The unit is one model, which is what a file of the folder holds. A model
is read with
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md),
passed to `f`, and written back to its own file through a temporary file
that replaces the old one once it is complete, so an interrupted step
leaves either the model as it was or the model as `f` made it, never
half of it. A file holding more than one model cannot be worked on this
way, since writing one of them back would have to rewrite the others,
and such a model is refused rather than risked.

The manifest is rebuilt from the specifications of the models that come
back, so a step that changes what a model is – its rank, its lag orders,
the chain it asks for – leaves the table describing the files. It
records what a model is, not how many draws are left in it, so thinning
the draws does not change it.

A step over a subset is a step over the whole folder applied to part of
it: the position a model has in the whole manifest is what is passed as
`index` and what numbers its seed, so a run taken model by model draws
what a run over all of them draws. That is what makes a long estimation
resumable, and the estimation steps for a folder pass `models` on.

With more than one core, `f` is sent to the workers as it stands, and
the workers load the installed package. A function that calls something
the package does not export, or something only a newer version of it
has, fails there rather than in the session it was written in. Pass what
such a function needs as an argument instead.

For the posterior itself, prefer
[`bayests_files`](https://franzmohr.github.io/bvartools/reference/bayests_files.md)
over a cluster: the draws are the largest thing a model has, and reading
them into R and writing them out again costs more than the sampler. An
idle R worker is not cheap either – on Windows a fresh one commits about
2 GB before it does anything and about 4 GB once it has loaded this
package.

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

folder <- file.path(tempdir(), "bvartools-example-map")
unlink(folder, recursive = TRUE)
dir.create(folder, recursive = TRUE)
write_to_hdf5(models, folder = folder)

stored <- open_models(folder)

# What add_priors() for a folder does
stored <- map_models(stored, function(model) {
  add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1))
})

# Reading only: how many observations each model was built on
map_models(stored, function(model) nrow(model$data$train$y), write = FALSE)
#> $`VarNormalWishart-varsel=none-p=01-n=1-001`
#> [1] 89
#> 
#> $`VarNormalWishart-varsel=none-p=02-n=1-001`
#> [1] 89
#> 
```
