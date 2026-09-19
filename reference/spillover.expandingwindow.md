# Rolling Spillover Index

Produces the total spillover index of Diebold and Yilmaz (2012) over the
windows of an object of class 'expandingwindow'.

## Usage

``` r
# S3 method for class 'expandingwindow'
spillover(object, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow' whose models carry posterior
  draws.

- ...:

  arguments passed forward to
  [`spillover.bvarmodel`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md).

## Value

A time-series object of class 'bvarspilloverts' with the median of the
total index and the bounds of its credible interval, one row per window.
Windows whose estimation failed are `NA`.

## Details

The index over a moving sample is the chart this literature is built
around, and
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
already produces the models it needs. Each window contributes one value,
dated at its own last observation, so the result reads as the
connectedness of the system as it was known up to that date.

Note that the windows expand rather than roll: every one of them starts
at the beginning of the sample, so later values are estimated on more
data. Diebold and Yilmaz use a fixed width instead, which trades that
away for a constant amount of information per point.

## References

Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive:
Predictive directional measurement of volatility spillovers.
*International Journal of Forecasting, 28*(1), 57–66.
[doi:10.1016/j.ijforecast.2011.02.006](https://doi.org/10.1016/j.ijforecast.2011.02.006)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100
e1 <- window(e1, end = c(1978, 4))

# Generate model data
model <- create_bvarmodel(e1, p = 1, deterministic = "const",
                          iterations = 50, burnin = 10)
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Estimate over expanding windows
windows <- use_expanding_window(model, start = c(1976, 1))
windows <- lapply(windows, function(w) {
  add_posterior_coefficients(add_initial_values(w))
})
class(windows) <- append("expandingwindow", class(windows))

# The index as the sample grows
si <- spillover(windows, n_ahead = 5)
plot(si)

```
