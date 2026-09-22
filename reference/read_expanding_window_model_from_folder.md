# Import Models from HDF5 Files

Imports model information and posterior draws from an HDF5 file.

## Usage

``` r
read_expanding_window_model_from_folder(folder)
```

## Arguments

- folder:

  Path to a folder with HDF5 files containing model data.

## Value

A list of class 'expandingwindow' with one model per file in `folder`,
each read by
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md),
in the order [`list.files`](https://rdrr.io/r/base/list.files.html)
returns the files.
