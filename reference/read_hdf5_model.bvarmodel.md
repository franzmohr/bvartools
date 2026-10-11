# Read a VAR or VEC Model from an Open File

The methods of
[`read_hdf5_model`](https://rdrr.io/pkg/bayests/man/read_hdf5_model.html)
that
[`read_model_from_hdf5`](https://rdrr.io/pkg/bayests/man/read_model_from_hdf5.html)
reaches for a file of a VAR or VEC model: one written by
[`write_to_hdf5`](https://rdrr.io/pkg/bayests/man/write_to_hdf5.html),
which records the class in `/model/rclass`, or by the BayesTS command
line, whose algorithm names it.

## Usage

``` r
# S3 method for class 'bvarmodel'
read_hdf5_model(x, h5_root, draws = NULL, ...)

# S3 method for class 'bvecmodel'
read_hdf5_model(x, h5_root, draws = NULL, ...)
```

## Arguments

- x:

  an empty object of class 'bvarmodel' or 'bvecmodel', the class to
  read.

- h5_root:

  the group of the open file the model hangs under.

- draws:

  the positions of the draws to read, or `NULL` for every draw.

- ...:

  not used.

## Value

An object of class 'bvarmodel' or 'bvecmodel' with the elements the file
holds: `model`, `data`, `priors`, `initial` and, where the model has
been estimated, `posterior`. The draws of a sampler are
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) objects; the per-period
posterior of a discounted model is a plain matrix, since its rows are
periods rather than draws.
[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the elements one by one.

## Details

The groups are read one by one rather than as a tree, since this layout
predates the tree: a time series is restored with its `tsp`, the scale
and centre of a rescaled error correction term are put back, and the
sign restrictions and dummy variables, which are tables, are read into
data frames.

## Examples

``` r

data("e1")
e1 <- diff(log(e1)) * 100
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

# The method behind read_model_from_hdf5() for VAR and VEC models
file <- file.path(tempdir(), "bvartools-example-read.h5")
unlink(file)
write_to_hdf5(model, filename = file)
stored <- read_model_from_hdf5(file)
class(stored)
#> [1] "bvarmodel" "list"     
```
