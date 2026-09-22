# Add Forecast Input Data

Generates and adds data matrices for forecast simulation to the elements
of an object of class 'expandingwindow'.

## Usage

``` r
# S3 method for class 'expandingwindow'
add_forecast_input(object, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow' containing objects that can be
  passed on to their respective `add_forecast_input` method.

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
model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

model <- use_expanding_window(model, start = 1982.25)

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add data used for forecast calculation
model <- add_forecast_input(model, n_ahead = 4)
```
