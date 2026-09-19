# Add Initial Values of an MCMC Chain

Adds initial values to a list of models by passing each element to the
respective method.

## Usage

``` r
# S3 method for class 'expandingwindow'
add_initial_values(object, ...)
```

## Arguments

- object:

  a list of class 'expandingwindow', usually, the output of a call to
  [`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md).

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with initial values added to each of its models,
as described in
[`add_initial_values.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md).

## Examples

``` r

data("us_macrodata")

# Starting period of the forecasting exercise
start_period <- 2007

# AR(1) models as benchmark
model <- create_bvarmodel(data = us_macrodata, p = 1,
                          deterministic = "none",
                          error = "gamma",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in draws should be much higher.

# Obtain objects for expanding window estimation
model <- use_expanding_window(model, start = start_period)

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(shape = 3, rate = 1))

# Add initial values
model <- add_initial_values(model)
```
