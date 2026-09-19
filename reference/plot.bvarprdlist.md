# Plotting Forecasts of a List of Bayesian VAR Models

A plot function for objects of class 'bvarprdlist'.

## Usage

``` r
# S3 method for class 'bvarprdlist'
plot(x, n_pre = NULL, ci = 0.95, ...)
```

## Arguments

- x:

  an object of class 'bvarprdlist', usually, a result of a call to
  [`predict.modellist`](https://franzmohr.github.io/bvartools/reference/predict.modellist.md).

- n_pre:

  number of plotted observations that precede the forecasts. If `NULL`
  (default), all available observations will be plotted.

- ci:

  interval used to calculate the credible bands of the forecasts.

- ...:

  further graphical parameters, which are passed on to
  [`plot.ts`](https://rdrr.io/r/stats/plot.ts.html).

## Details

The forecasts of the models in `x` are plotted on a grid, where each
column corresponds to a model and each row to an endogenous variable.
