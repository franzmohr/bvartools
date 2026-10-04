# Turn the Tree of a Model File into a Model

The step of
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
that turns what
[`read_bayests_tree`](https://franzmohr.github.io/bvartools/reference/write_bayests_tree.md)
reads into the object of a model package.

## Usage

``` r
from_bayests_tree(tree, ...)

# Default S3 method
from_bayests_tree(tree, ...)
```

## Arguments

- tree:

  the tree of a model file, as
  [`read_bayests_tree`](https://franzmohr.github.io/bvartools/reference/write_bayests_tree.md)
  returns it, with the class recorded in the file's `/model/rclass`.

- ...:

  further arguments passed to or from other methods.

## Value

An object of the class of `tree`, as its package defines it.

## Details

[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
reads the VAR and VEC models of this package itself. For a file whose
`rclass` names any other class it reads the tree, gives it that class
and calls this generic, so a package that writes its models with
[`write_bayests_tree`](https://franzmohr.github.io/bvartools/reference/write_bayests_tree.md)
reads them back by registering a method for its class. The method undoes
the translation its writer made: names, orderings and shapes that
BayesTS wants and the package's object does not.

The package named in the file's `/model/rpackage` is loaded first, if it
is installed, so that its methods are registered whether or not it is
attached. A file of the BayesTS command line carries no such attribute;
one of a factor model is read with dfmtools loaded.

Without a method the default refuses, and names the class it found no
method for.

**A method keeps the attributes of `/model` it does not know in element
`model` of what it returns.**
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
records there where a model stood in the 'modellist' it was written
from, and
[`read_models_from_folder`](https://franzmohr.github.io/bvartools/reference/read_models_from_folder.md)
rebuilds the list from it.
