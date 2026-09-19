# Combine Models

Convenience function that combines multiple models into one object.

## Usage

``` r
combine_models(...)
```

## Arguments

- ...:

  One or multiple lists that contain either model objects or posterior
  objects.

## Value

A list of class 'modellist' that holds the models of all arguments in
one flat list, in the order in which they were given.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

# Load data
data(e6)

# Create multiple VAR models
model_1 <- create_bvarmodel(diff(e6), p = 0:2, deterministic = "const")

# Create multiple VEC models
model_2 <- create_bvecmodel(e6, p = 1, r = 0:1, const = "unrestricted")

# Combine the models in one object
model <- combine_models(model_1, model_2)
```
