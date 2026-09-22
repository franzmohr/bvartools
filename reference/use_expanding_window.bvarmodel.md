# Expanding Window Estimation

Creates objects for expanding window posterior simulation.

## Usage

``` r
# S3 method for class 'bvarmodel'
use_expanding_window(object, start, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel' containing model specifications and
  input data. Usually, the output of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- start:

  the start period of the prediction of the first iteration of the
  expanding window approach.

- ...:

  arguments passed forward to method.

## Value

A list of class 'expandingwindow' with one object of class 'bvarmodel'
per window. The training sample of the first ends in the period before
`start`, and each further window adds one period. Posterior draws,
starting values and forecast input that `object` already carries belong
to the whole sample and are not copied into the windows; a warning says
so.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

data("us_macrodata")

# Starting period of the forecasting exercise
start_period <- 2007

# Create model
model <- create_bvarmodel(data = us_macrodata,
                          p = 1,
                          deterministic = "none",
                          seasonal = FALSE,
                          tvp = FALSE,
                          error = "gamma",
                          iterations = 10,
                          burnin = 2)
# Chosen number of iterations and burn-in draws should be much higher.

# Create multiple model objects for expanding window
model <- use_expanding_window(model, start = start_period)
```
