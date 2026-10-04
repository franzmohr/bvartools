# Sign and Zero Restricted Rotations of Posterior Draws

The identification step of
[`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md)
on its own: rotations of the reduced-form draws of a VAR that satisfy
sign and zero restrictions on the impulse responses, with the importance
weights of Arias, Rubio-Ramirez and Waggoner (2018). It is exported for
packages that estimate a VAR as part of something larger – the state
equation of a factor augmented VAR, for example – and want to identify
it the same way.

## Usage

``` r
arias_rubio_ramirez_waggoner_2018(
  draws,
  restrictions,
  variables,
  lags,
  max_tries = 1,
  smooth = TRUE,
  one_sided = FALSE
)
```

## Arguments

- draws:

  a list with one element per posterior draw, each a list with elements
  `A`, the \\K \times M\\ coefficient matrix with the \\p\\ lag blocks
  first and any deterministic terms after them, and `Sigma`, the \\K
  \times K\\ error covariance.

- restrictions:

  a data frame of the restrictions, as in
  [`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md):
  columns `impulse`, `response`, `sign` and optionally `horizon`, the
  first two naming elements of `variables`. A table without a zero
  restriction is accepted here; the rotations are then drawn uniformly
  and the weights are equal.

- variables:

  the names of the \\K\\ variables, in the order of the rows of `A`.
  Shocks are named after them and their columns are drawn in this order.

- lags:

  the lag order \\p\\.

- max_tries:

  integer. The largest number of rotations drawn for each posterior
  draw. See 'Details' of
  [`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md).

- smooth:

  logical. Should the importance weights be Pareto smoothed?

- one_sided:

  logical. Should the numerical derivative behind the importance weights
  be taken on one side only?

## Value

A list with

- q:

  a matrix with one row per draw holding its rotation, by column, and
  `NA` for a draw with no admissible rotation;

- weights:

  the normalised importance weights, zero for such a draw;

- restrictions:

  the restrictions, with positions in place of names;

- tries:

  the number of rotations drawn in all;

- accepted:

  the number of draws with an admissible rotation;

- effective_sample_size, max_weight_share, pareto_k:

  the diagnostics of the importance sample, as
  [`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md)
  records them.

## Details

The rotation \\Q\\ of a draw maps into the impact matrix \\L Q\\, with
\\L\\ the lower Cholesky factor of `Sigma`: the response of the
variables to the shocks on impact.

The function draws from the random number generator in the order
[`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md)
does, so the same seed gives the same rotations through either. It warns
when the importance sample is not fit to summarise, and stops when no
draw has an admissible rotation.

## References

Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference
based on structural vector autoregressions identified with sign and zero
restrictions: Theory and applications. *Econometrica, 86*(2), 685-720.
[doi:10.3982/ECTA14468](https://doi.org/10.3982/ECTA14468)

## See also

[`add_sign_zero_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_zero_restrictions.md),
which applies it to a 'bvarmodel' and resamples the draws.
