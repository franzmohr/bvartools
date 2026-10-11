# Time Series Windows

Restricts the observations in the training sample of a model of class
'bvarmodel' to the specified time window.

## Usage

``` r
# S3 method for class 'bvarmodel'
window(x, start = NULL, end = NULL, ...)
```

## Arguments

- x:

  an object of class 'bvarmodel'.

- start:

  the start time of the period of interest.

- end:

  the end time of the period of interest.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'bvarmodel'.

## Details

Posterior draws that form a path over the periods of the training sample
are cut to the periods that remain, so that they still refer to the same
observations as the data. These are the coefficients of a model with
time varying parameters, the draws of the error term that vary by period
– as under stochastic volatility – the panel completed by a model with
`missing = "estimate"`, and the pointwise log-likelihood. Draws that do
not vary by period are left as they are. The posterior of a discounted
model, which holds the moments of every period rather than draws, is cut
to the same periods. In either case the draws or moments that remain
were estimated from the whole sample; re-estimate the model to condition
on the window alone.

The period of a volatility break moves with the window. A window that
ends before the break drops it, with its prior and its draws, and a
warning; one that starts after it is refused, since a break has to fall
inside the sample.

**A window that ends before the sample drops the forecast input, the
forecasts and the realised values they were scored against**, with a
warning: they start after the original end, from its last observations.
Call
[`add_forecast_input()`](https://rdrr.io/pkg/bayests/man/add_forecast_input.html)
again to forecast from the end of the window.

## Examples

``` r

data("e1")
e1 <- diff(log(e1)) * 100
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

# The last ten years of the sample
recent <- window(model, start = c(1972, 1))
start(recent$data$train$y)
#> [1] 1972    1
```
