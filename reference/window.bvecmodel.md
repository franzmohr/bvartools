# Time Series Windows

Restricts the observations in the training sample of a model of class
'bvecmodel' to the specified time window.

## Usage

``` r
# S3 method for class 'bvecmodel'
window(x, start = NULL, end = NULL, ...)
```

## Arguments

- x:

  an object of class 'bvecmodel'.

- start:

  the start time of the period of interest.

- end:

  the end time of the period of interest.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'bvecmodel'.

## Details

Posterior draws that form a path over the periods of the training sample
are cut to the periods that remain, so that they still refer to the same
observations as the data. These are the coefficients and cointegration
vectors of a model with time varying parameters, the draws of the error
term that vary by period – as under stochastic volatility – and the
pointwise log-likelihood. Draws that do not vary by period, such as
those of \\\rho\\, are left as they are.

## Examples

``` r

data("e6")
model <- create_bvecmodel(e6 * 100, p = 1, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burn-in should be much higher.
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

recent <- window(model, start = c(1980, 1))
start(recent$data$train$y)
#> [1] 1980    1
```
