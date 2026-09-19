# Export to HDF5 File

A generic function used to export the content of models into HDF5 files.
The function invokes particular methods which depend on the class of the
first argument.

## Usage

``` r
write_to_hdf5(object, ...)
```

## Arguments

- object:

  an object of a class, for which a method should be called.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`write_to_hdf5.bvarmodel`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.bvarmodel.md),
[`write_to_hdf5.bvecmodel`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.bvecmodel.md),
[`write_to_hdf5.expandingwindow`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.expandingwindow.md),
[`write_to_hdf5.modellist`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.modellist.md).
