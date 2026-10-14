# The prior of Student t errors

Checks `sigma$student_t` against the model and returns the prior the
degrees of freedom are drawn under, or `NULL` for a model with normal
errors.

## Usage

``` r
student_t_prior(object, sigma)
```

## Arguments

- object:

  a model object whose `model$student_t` says whether its errors are
  Student t.

- sigma:

  the list given to
  [`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html), whose
  element `student_t` is either `list(df = )`, holding the degrees of
  freedom fixed, or `list(shape = , rate = )`, the gamma prior of Juarez
  and Steel (2010) under which they are drawn.

## Value

A named list of the prior, or `NULL`.

## Details

Exported so that a package building its own models checks the prior the
way [`add_priors`](https://rdrr.io/pkg/bayests/man/add_priors.html)
does. Like every other hyperparameter it has no default, and a prior
given to a model with normal errors is refused rather than ignored.

## References

Juarez, M. A. and Steel, M. F. J. (2010). Model-based clustering of
non-Gaussian panel data based on skew-t distributions. *Journal of
Business & Economic Statistics*, 28(1), 52–66.
