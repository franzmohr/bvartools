# Transform a VEC Model to a VAR in Levels

A generic function used to transform a vector error correction model
into its VAR form. The function invokes particular methods which depend
on the class of the first argument.

## Usage

``` r
vec_to_var(object, ...)
```

## Arguments

- object:

  an object with suitable input data passed forward to method.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`vec_to_var.bvecmodel`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md),
[`vec_to_var.modellist`](https://franzmohr.github.io/bvartools/reference/vec_to_var.modellist.md).
