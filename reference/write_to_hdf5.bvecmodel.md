# Export to HDF5 File

Exports the content of an object of class 'bvarmodel' to an HDF5 file.

## Usage

``` r
# S3 method for class 'bvecmodel'
write_to_hdf5(object, filename, group = "", ...)
```

## Arguments

- object:

  list of class 'bvarmodel'.

- filename:

  path to the file, in which output should be stored.

- group:

  the group the model's tree should hang under inside its file. Defaults
  to `""`, the root of the file, which is where a file holding a single
  model puts it. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The path to the written file, invisibly.

## Details

With a `group` every path is written under it instead of at the root, so
one file can hold several models side by side. The spelling is the one
the BayesTS command line uses for its `--group` flag: a leading slash
and no trailing slash, with `""` for the root. Intermediate groups are
created as needed, and
[`list_models_in_hdf5`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md)
reports which groups of a file hold a model.

What must not already be there is the model rather than the file.
Without a `group` the model is the whole file, so an existing file is
refused; with one, only that group has to be free, which is what lets a
second model be added beside the first.

A write that cannot be completed raises the error rather than absorbing
it, and undoes what the call created: a file it made is removed, a group
it added to an existing file is unlinked on its own so the models beside
it survive. A retry then meets the original problem rather than the
leftovers.

## Examples

``` r

# Load data 
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Save model
path_to_model <- tempfile(fileext = ".h5")
write_to_hdf5(model, filename = path_to_model)
```
