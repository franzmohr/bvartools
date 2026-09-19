# Scale Error Correction

Scales, centres, or centres and scales the series in the error
correction term of an object of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
scale_error_correction(object, scale = TRUE, centre = FALSE, ...)
```

## Arguments

- object:

  object of class 'bvecmodel'.

- scale:

  logical. If `TRUE`, the default, the series are divided by scaling
  factors. See 'Details'.

- centre:

  logical. If `TRUE`, the sample mean of each stochastic series is
  subtracted from it. Default is `FALSE`. See 'Details'.

- ...:

  arguments passed forward to method.

## Value

An object of class 'bvecmodel'.

## Details

The function transforms element `object$data$train$w`, and
[`rescale_error_correction`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
undoes the transformation once the posterior draws are in.

With `scale = TRUE`, stochastic series are divided by the standard
deviation of the corresponding differenced series. If the time-series
object contains a column named `"trend"`, this series is divided by its
own standard deviation, i.e. in levels. The scaling factors are stored
as a new attribute of `object$data$train$w` named `"scale"`.

With `centre = TRUE`, the sample mean of each stochastic series is
subtracted from it. The deterministic terms restricted to the
cointegration space, which are the last columns of the term, are left as
they are. The means are stored as attribute `"centre"`, with zeros for
the deterministic terms. If both arguments are `TRUE`, the series are
centred first and scaled second, so that \\w_t\\ becomes \\D^{-1} (w_t -
m)\\.

Centring requires an unrestricted constant, which takes up what the
series lose: \\\Pi w_t + c = \Pi (w_t - m) + (c + \Pi m)\\. For a model
with constant coefficients it is therefore a reparameterisation that
leaves the likelihood and the priors on \\\alpha\\ and \\\beta\\ as they
are and changes only what the prior of the constant refers to. For a
model with time varying cointegration vectors it changes more. A step
\\\eta_t\\ of their state equation moves \\\beta_t^{\prime} w_t\\ by
\\\eta_t^{\prime} w_t\\, which for series far from zero acts as a random
walk intercept that no prior on the deterministic terms controls; see
section 'Prior on the cointegration space' of
[`cointspace_prior`](https://franzmohr.github.io/bvartools/reference/cointspace_prior.md).
On centred series the same step moves the term by \\\eta_t^{\prime}
(w_t - m)\\, and the intercept is left to the constant and its own
prior.

Starting values of the constant in `object$initial` are shifted by \\\Pi
m\\, computed from the starting values of the loadings and of \\\beta\\,
so that the chain starts from the model it simulates.

Neither transformation can be applied twice, and neither can be added to
a model that already carries the other: `rescale_error_correction` has
to be called in between.
