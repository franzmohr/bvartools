# Dynamic Multipliers of a VEC Model with Exogenous Variables

Computes the response of an endogenous variable to a change in a weakly
exogenous variable of an object of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
multipliers(x, ...)
```

## Arguments

- x:

  an object of class 'bvecmodel' with at least one exogenous variable.

- ...:

  arguments of
  [`multipliers.bvarmodel`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
  which does the work.

## Value

A time-series object of class 'bvarirf', which is what
[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md) returns,
so that the same `plot` method applies.

## Details

A multiplier is a statement about the levels of the endogenous
variables, so the model is put into its levels form with
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
and the multipliers of that form are returned. The transformation is
exact, and it is where the error correction term does its work: the
long-run relations enter the levels coefficients, so a change in an
exogenous variable that belongs to a cointegrating relation moves the
endogenous variables permanently, while at rank zero it moves them only
through the short-run terms.

The responses are therefore in the units of the levels of the endogenous
variables, whatever the error correction form has on its left-hand side.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md),
[`multipliers.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

data("e6")
set.seed(1)
exogen <- ts(cbind(g = as.numeric(e6[, "R"]) * 0.3 + rnorm(nrow(e6), sd = 0.1)),
             start = start(e6), frequency = frequency(e6))

model <- create_bvecmodel(data = e6, exogen = exogen, p = 2, s = 1, r = 1,
                          const = "unrestricted",
                          iterations = 100, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = 3, scale = 0.0001))
model <- add_posterior_coefficients(add_initial_values(model))

dm <- multipliers(model, impulse = "g", response = "R", n_ahead = 8)
```
