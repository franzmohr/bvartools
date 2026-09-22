# Add Forecast Errors

Generic function used to calculate forecast errors and add them to a
model object.

## Usage

``` r
add_forecast_errors(object, test_sample = NULL, ...)
```

## Arguments

- object:

  an object of a class, for which a method should be called.

- test_sample:

  a time-series object used as test data. If `NULL` (default), the
  values in `data$test$y` of the object are used, which is what a model
  carries after it has been scored once and what a model read from a
  file was written with.

- ...:

  arguments passed forward to method.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## See also

Methods:
[`add_forecast_errors.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.expandingwindow`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.expandingwindow.md),
[`add_forecast_errors.modellist`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.modellist.md).
