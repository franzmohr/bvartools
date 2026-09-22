# Plotting Forecasts of BVAR Models

A plot function for objects of class 'bvarprd'.

## Usage

``` r
# S3 method for class 'bvarprd'
plot(x, n_pre = NULL, ci = 0.95, ...)
```

## Arguments

- x:

  an object of class 'bvarprd', usually, a result of a call to
  [`predict.bvarmodel`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md).

- n_pre:

  number of plotted observations that precede the forecasts. If `NULL`
  (default), all available observations will be plotted.

- ci:

  interval used to calculate the credible bands of the forecasts.

- ...:

  further graphical parameters. Arguments `main`, `ylab`, `lty` and
  `plot.type` overwrite the defaults of the function.

## Value

`x`, invisibly. The function is called for its side effect, the plot.

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

# Plot forecasts
plot(pred)



```
