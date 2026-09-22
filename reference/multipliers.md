# Dynamic Multipliers

A generic function used to calculate the dynamic multipliers of a model
with weakly exogenous variables. The function invokes particular methods
which depend on the class of the first argument.

## Usage

``` r
multipliers(x, ...)
```

## Arguments

- x:

  an object with suitable input data passed forward to method.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `x`, as described on
the pages of the methods.

## Details

A dynamic multiplier is the response of the endogenous variables of a
model to a change in one of its weakly exogenous variables. It is what
an impulse response is for a shock to an error term, for a variable the
model does not explain: the exogenous variable is moved by hand and the
endogenous ones are followed.

The quantity is central to models whose foreign block is weakly
exogenous, such as the country models of a global VAR (Pesaran,
Schuermann and Weiner, 2004), where the responses to a change in the
foreign variables are what the country model has to say about the rest
of the world without the rest of the world being solved.

## References

Pesaran, M. H., Schuermann, T., Weiner, S. M. (2004). Modeling regional
interdependencies using a global error-correcting macroeconometric
model. *Journal of Business & Economic Statistics, 22*(2), 129-162.

## See also

Methods:
[`multipliers.bvarmodel`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md).

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`multipliers.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)
