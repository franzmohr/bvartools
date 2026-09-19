# Expanding Window Estimation

Creates objects for expanding window posterior simulation.

## Usage

``` r
# S3 method for class 'modellist'
use_expanding_window(object, start, ...)
```

## Arguments

- object:

  a list of objects containing model specifications and input data.
  Usually, the output of a call to a model creation function such as
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- start:

  the start period of the prediction of the first iteration of the
  expanding window approach.

- ...:

  arguments passed forward to method.

## Value

A list of class 'modellist', which consists of objects of class
'expandingwindow'.

## Examples

``` r

data("us_macrodata")

# Starting period of the forecasting exercise
start_period <- 2007

# Create model
model <- create_bvarmodel(data = us_macrodata,
                          p = 1:4,
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
