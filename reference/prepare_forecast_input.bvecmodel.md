# Prepare Forecast Input

Prepares the regressors of the forecast periods of a VEC model, in
levels.

## Usage

``` r
# S3 method for class 'bvecmodel'
prepare_forecast_input(object, ...)
```

## Arguments

- object:

  an object of class 'bvecmodel'.

- ...:

  arguments passed to
  [`prepare_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.bvarmodel.md),
  such as `n_ahead`, `deterministic` and `exogen`.

## Value

A list with the forecast horizon in element `h` and the regressors of
the forecast periods in element `x`, as returned by
[`prepare_forecast_input.bvarmodel`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.bvarmodel.md).

## Details

A VEC model is forecast in levels, so its forecast regressors are those
of its VAR representation. They are prepared from the data and
specification
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
builds; the posterior draws are not converted.

The forecast of a VEC model starts right after its estimation sample,
and argument `y_start` of the 'bvarmodel' method is refused: the VEC
samplers forecast from their own draws and the end of the sample, and a
time varying model steps its states from there. To start a forecast from
values after the sample, convert the model with
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
and give `y_start` to the resulting 'bvarmodel'.
