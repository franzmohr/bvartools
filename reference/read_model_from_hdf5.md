# Import Models from HDF5 Files

Imports model information and posterior draws from an HDF5 file.

## Usage

``` r
read_model_from_hdf5(filename, group = "")
```

## Arguments

- filename:

  Path to an HDF5 file containing model data.

- group:

  the group the model's tree hangs under inside its file. Defaults to
  `""`, the root of the file, which is where a file holding a single
  model puts it. See 'Details'.

## Details

With a `group` every path the reader looks for is read under it instead
of at the root, so one file can hold several models side by side.
[`list_models_in_hdf5`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md)
reports which groups of a file hold one.

The spelling is the one the BayesTS command line uses for its `--group`
flag: a leading slash and no trailing slash, with `""` for the root.
