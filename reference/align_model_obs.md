# Align Observations Across Models

Generic function used to restrict each model in a list to the set of
observations common to all models, ensuring that comparisons are
computed on the same underlying sample.

## Usage

``` r
align_model_obs(object, ...)
```

## Arguments

- object:

  a list of models.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`align_model_obs.modellist`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md).
