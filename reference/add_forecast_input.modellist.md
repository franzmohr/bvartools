# Add Forecast Input Data

Generates and adds data matrices for forecast simulation to the elements
of an object of class 'modellist'.

## Usage

``` r
# S3 method for class 'modellist'
add_forecast_input(object, ...)
```

## Arguments

- object:

  an object of class 'modellist' containing objects that can be forward
  to their respective \`add_forecast_input\` method.

- ...:

  arguments passed forward to method.

## Value

The object in `object` with forecast input added to each of its models,
as described in
[`add_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md).

## Examples

``` r

# Load data
data("e1")
orig <- diff(log(e1)) * 100
train <- window(orig, end = c(1982, 2))


# Create model
model <- create_bvarmodel(data = train, p = 0:2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

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
```
