# Import Models from a Folder of HDF5 Files

Imports every model stored below a folder.

## Usage

``` r
read_models_from_folder(folder, draws = NULL)
```

## Arguments

- folder:

  Path to a folder containing model data.

- draws:

  the draws to read of every model, as
  [`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
  takes them: `NULL`, the default, for the whole chain, a vector of
  positions for those draws, or `integer(0)` for none of them.

## Value

A named list of class 'modellist', or of class 'expandingwindow' if the
models say they belong to one.

## Details

The folder is walked recursively, and each HDF5 file in it is asked
which of its groups hold a model – so a file holding several models
contributes all of them. See
[`list_models_in_hdf5`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md)
for what counts as one.

The result is a flat named list rather than a nesting that mirrors the
directory tree. Each name is the path of the file relative to `folder`
without its extension, followed by the group where a model does not sit
at the root of its file. Names are what a caller needs to tell one model
from another – which sub-model of a global model it is, say – and a
nesting whose depth depended on where the caller pointed could not
provide them.

The one exception is an expanding window: the windows that
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
wrote to one directory come back as one element of class
'expandingwindow', named after that directory, so that the windows of
different specifications are not pooled. A folder that holds a single
expanding window is returned as that expanding window. A model list
written by
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
comes back in the order it was written in, which the files record; other
models are in the order of their names.

The class of each model comes from the `rclass` attribute the writer
records, not from its file name.

`draws` is what makes a folder larger than the session readable: with
`integer(0)` the models come back with their specification, their data
and their priors and no draws, and with a vector of positions the chain
is read a piece at a time. For working on the models rather than reading
them,
[`open_models`](https://franzmohr.github.io/bvartools/reference/open_models.md)
hands them over one at a time instead.

## Examples

``` r

# Load data
data("e1")
train <- diff(log(e1)) * 100

# Create and store models
folder <- file.path(tempdir(), "models")
dir.create(folder, showWarnings = FALSE)
for (p in 1:2) {
  model <- create_bvarmodel(data = train, p = p, deterministic = "const",
                            iterations = 10, burnin = 10)
  write_to_hdf5(model, filename = file.path(folder, paste0("model-", p, ".h5")))
}

models <- read_models_from_folder(folder)
names(models)
#> [1] "model-1" "model-2"
```
