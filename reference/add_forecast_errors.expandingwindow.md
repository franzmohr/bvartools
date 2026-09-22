# Add Forecast Errors

Calculates and adds forecast errors for a list of Bayesian models.

## Usage

``` r
# S3 method for class 'expandingwindow'
add_forecast_errors(object, test_sample = NULL, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow'.

- test_sample:

  a time-series object used as test data. If `NULL` (default), the
  values in `data$test$y` of the object are used.

- ...:

  arguments passed forward to method.

## Value

The object in `object` with forecast errors added to each of its models,
as described in
[`add_forecast_errors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md).

## Examples

``` r

# Load data
data("e1")
orig <- diff(log(e1)) * 100
train <- window(orig, end = c(1982, 2))


# Create model
model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

model <- use_expanding_window(model, start = 1982.25)

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Add data used for forecast calculation
model <- add_forecast_input(model, n_ahead = 4)

# Add forecasts
model <- add_posterior_forecasts(model)

# Add forecast errors
model <- add_forecast_errors(model, test_sample = orig)
```
