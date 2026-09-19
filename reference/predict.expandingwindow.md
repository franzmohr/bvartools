# Predict Method for Objects of Class expandingwindow

Forecasting Bayesian VAR objects in a list of class 'expandingwindow'.

## Usage

``` r
# S3 method for class 'expandingwindow'
predict(object, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow'.

- ...:

  additional arguments.

## Value

A time-series object of class 'expandwindbvarprdlist'.

## Examples

``` r

data("us_macrodata")

model <- create_bvarmodel(data = us_macrodata,
                          p = 1,
                          deterministic = "none",
                          error = "gamma",
                          iterations = 10,
                          burnin = 2)
# Chosen number of iterations and burn-in draws should be much higher.

# Obtain objects for expanding window estimation
model <- use_expanding_window(model, start = 2007)

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(shape = 3, rate = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

model <- add_forecast_input(model, n_ahead = 10)
model <- add_posterior_forecasts(model)

# Predict
prd <- predict(model, n_ahead = 4)
```
