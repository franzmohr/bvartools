# Spillover Index

Guard method for objects of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
spillover(object, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel'.

- ...:

  additional arguments.

## Value

Nothing. The method raises an error, since spillover indices of a VEC
model are obtained from its VAR representation: apply
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
first and then
[`spillover`](https://franzmohr.github.io/bvartools/reference/spillover.md)
to the resulting 'bvarmodel'.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`add_sign_zero_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`historical_decomposition()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.md),
[`historical_decomposition.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/historical_decomposition.bvarmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`split_quantile_grid()`](https://franzmohr.github.io/bvartools/reference/split_quantile_grid.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)
