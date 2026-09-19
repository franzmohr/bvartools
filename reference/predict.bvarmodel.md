# Predict Method for Objects of Class bvar

Forecasting a Bayesian VAR object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
predict(object, n_ahead = NULL, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel'.

- n_ahead:

  number of steps ahead at which to predict. If `NULL` (default), every
  period simulated by
  [`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md),
  i.e. the horizon given to
  [`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md).

- ...:

  additional arguments.

## Value

An object of class 'bvarprd', a list with element `fcst`, an \\h \times
K \times S\\ array of the \\S\\ forecast draws, whose first two
dimensions are named after the forecast periods and the endogenous
variables, and element `y`, the endogenous variables of the training
sample.

## Details

For the VAR model \$\$A_0 y_t = \sum\_{i = 1}^{p} A\_{i} y\_{t-i} +
\sum\_{i = 0}^{s} B\_{i} x\_{t-i} + C D_t + u_t,\$\$ with \\u_t \sim
N(0, \Sigma)\\ the function produces `n_ahead` forecasts.

## References

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100
e1 <- window(e1, end = c(1978, 4))

# Generate model data
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 10)
# Chosen number of iterations and burnin should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# Add the data the forecasts are produced from
model <- add_forecast_input(model, n_ahead = 10)

# Simulate the forecasts
model <- add_posterior_forecasts(model)

# Collect them
pred <- predict(model)
```
