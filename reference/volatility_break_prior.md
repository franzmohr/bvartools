# The prior of a volatility break

Checks `sigma$volatility` against the model and returns the prior a
volatility break is drawn under, or `NULL` for a model without one.

## Usage

``` r
volatility_break_prior(object, sigma)
```

## Arguments

- object:

  a model object whose `model$volatility_break` says whether it has a
  break.

- sigma:

  the list given to
  [`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html), whose
  element `volatility` is the prior: `shape`, the Pareto shape of every
  free scale, and `rho_shape1` and `rho_shape2`, the two shapes of the
  beta prior on the rate the excess decays at.

## Value

A named list of the prior, or `NULL`.

## Details

Exported for the same reason as
[`check_volatility_break`](https://franzmohr.github.io/bvartools/reference/check_volatility_break.md):
a package that builds its own models checks the prior the same way.
There is no default, as there is none for the rest of `sigma`, and the
message names the values of Lenza and Primiceri (2022). A prior given to
a model without a break is refused rather than ignored, and a break
without a prior is refused as well.

## See also

[`check_volatility_break`](https://franzmohr.github.io/bvartools/reference/check_volatility_break.md)
