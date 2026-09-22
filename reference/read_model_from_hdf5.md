# Import Models from HDF5 Files

Imports model information and posterior draws from an HDF5 file.

## Usage

``` r
read_model_from_hdf5(filename, group = "", draws = NULL)
```

## Arguments

- filename:

  Path to an HDF5 file containing model data.

- group:

  the group the model's tree hangs under inside its file. Defaults to
  `""`, the root of the file, which is where a file holding a single
  model puts it. See 'Details'.

- draws:

  the draws to read, as their positions in the chain. Defaults to
  `NULL`, every draw. An integer vector reads those draws of every block
  of the posterior, and `integer(0)` reads none of them, which gives the
  model, its data and its priors without the draws. See 'Details'.

## Value

An object of class 'bvarmodel' or 'bvecmodel', depending on the
algorithm recorded in the file, with the elements that file holds:
`model`, `data`, `priors`, `initial` and, where the model has been
estimated, `posterior`. The draws of a sampler are
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) objects; the per-period
posterior of a discounted model is a plain matrix, since its rows are
periods rather than draws.
[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the elements one by one.

## Details

With a `group` every path the reader looks for is read under it instead
of at the root, so one file can hold several models side by side.
[`list_models_in_hdf5`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md)
reports which groups of a file hold one.

The spelling is the one the BayesTS command line uses for its `--group`
flag: a leading slash and no trailing slash, with `""` for the root.

`draws` reads part of a chain. Every block of a posterior holds one row
per draw, and only the rows asked for are read from the file, so a
caller that works through a long chain in pieces – solving a global
model draw by draw, for instance – never holds more of it than the piece
it is working on. `integer(0)` reads a model without its draws, which is
what a step needs that only looks at what a model is.

A partial read cannot describe the chain it came from, so the blocks it
returns are labelled as a chain of their own, from one to the number of
draws read. A full read keeps the labels the file carries.
