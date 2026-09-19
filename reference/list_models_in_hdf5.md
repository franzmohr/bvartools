# Models in an HDF5 File

Finds the groups of an HDF5 file that hold a model.

## Usage

``` r
list_models_in_hdf5(filename, group = "")
```

## Arguments

- filename:

  path to an HDF5 file.

- group:

  the group to search under. Defaults to `""`, the root of the file,
  which searches all of it.

## Value

A character vector of group names, sorted, empty if the file holds no
model. A model at the root of its file is reported as `""`, which is the
value
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
takes for it.

## Details

A model is a group with a `model` subgroup carrying an `algorithm`
attribute. The search stops at a model rather than descending into its
`data`, `priors` and `posterior`, which are its own subtree and not a
place further models could be.

This is the same rule the BayesTS command line applies for its
`--all-groups` flag, so both agree on which groups of a file are models.

## Examples

``` r

# Load data
data("e1")
train <- diff(log(e1)) * 100

# Create and store two models in one file
path_to_model <- tempfile(fileext = ".h5")
for (p in 1:2) {
  model <- create_bvarmodel(data = train, p = p, deterministic = "const",
                            iterations = 10, burnin = 10)
  write_to_hdf5(model, filename = path_to_model,
                group = paste0("/models/", p))
}

list_models_in_hdf5(path_to_model)
#> [1] "/models/1" "/models/2"
```
