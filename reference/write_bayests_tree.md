# Write and Read the Tree of a BayesTS Model File

Writes a nested list into an HDF5 file as a model's tree, or reads such
a tree back as a nested list. These are the steps beneath
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
and
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md),
exported for packages that keep models of their own in BayesTS files –
the dynamic factor models of dfmtools, for instance.

## Usage

``` r
write_bayests_tree(tree, filename, group = "")

read_bayests_tree(filename, group = "", draws = NULL)
```

## Arguments

- tree:

  a named list. An element that is itself a list is a group; any other
  element is a dataset. See 'Details'.

- filename:

  path to the HDF5 file.

- group:

  the group the tree hangs under inside its file. Defaults to `""`, the
  root of the file, which is where a file holding a single model puts
  it. The spelling is that of the BayesTS command line's `--group` flag:
  a leading slash and no trailing slash.

- draws:

  the draws to read from the blocks of `/posterior`, as their positions
  in the chain. Defaults to `NULL`, every draw; see
  [`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md).

## Value

`write_bayests_tree()` returns the path to the written file, invisibly.
`read_bayests_tree()` returns the tree as a nested list.

## Details

A model package translates its object into the layout BayesTS reads –
`/model`, `/data`, `/priors`, `/initial` and `/posterior`, named and
shaped as the sampler expects – and hands the result to
`write_bayests_tree()`. What is left to this function is what every
model shares: refusing to overwrite a model, undoing a write that fails
half-way, and the attributes that let a value be read back as what it
was. **It does not check the tree against any sampler**; a file that
BayesTS refuses is the translation's mistake, and `bayests check` is the
way to find it.

The tree is mapped onto the file as follows:

- An element that is a list becomes a group of the same name, and its
  elements go below it. `NULL` elements are left out.

- The element `.attributes` of a list is not a group but the attributes
  of the group it is in: a named list of values, each written as an
  attribute. This is where `/model` keeps the specification –
  `algorithm`, `k`, `p` and the rest – and where `rclass`, the class
  [`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
  returns, belongs, and `rpackage`, the package that defines that class,
  whose namespace the reader loads for its methods.

- Any other element becomes a dataset. HDF5 stores dimensions in the
  reverse of R's order, so a matrix of `tt` rows and `k` columns is a
  `(k, tt)` dataset, which is what BayesTS calls `(tt, k)` "on paper". A
  vector is a one-dimensional dataset and is marked so that it is read
  back as a vector. A time series carries its column names and `tsp`,
  and draws of class [`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html)
  their start, end and thinning interval, so that both are read back as
  what they were.

Without a `group` the tree is the whole file, so an existing file is
refused; with one, only that group has to be free, which is how a second
model is added beside the first. A write that cannot be completed raises
the error and undoes what the call created.

`read_bayests_tree()` is the inverse, for any model file, including one
written by the BayesTS command line: groups become lists, the attributes
of each group its element `.attributes`, and datasets values, restored
through the attributes above where the file carries them and as matrices
where it does not.

A file whose `/model` carries `rclass` is read by
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
through
[`from_bayests_tree`](https://franzmohr.github.io/bvartools/reference/from_bayests_tree.md),
which is the method a model package provides to turn the tree back into
its object.

## See also

[`from_bayests_tree`](https://franzmohr.github.io/bvartools/reference/from_bayests_tree.md),
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md),
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md).

## Examples

``` r

path <- tempfile(fileext = ".h5")
tree <- list(
  "model" = list(".attributes" = list("algorithm" = "Example", "k" = 2L)),
  "data" = list("train" = list("y" = matrix(rnorm(20), 10, 2)))
)
write_bayests_tree(tree, filename = path)
str(read_bayests_tree(path))
#> List of 2
#>  $ data :List of 1
#>   ..$ train:List of 1
#>   .. ..$ y: num [1:10, 1:2] 1.2609 -1.0503 -0.3901 0.0965 0.3258 ...
#>  $ model:List of 1
#>   ..$ .attributes:List of 2
#>   .. ..$ algorithm: chr "Example"
#>   .. ..$ k        : int 2
```
