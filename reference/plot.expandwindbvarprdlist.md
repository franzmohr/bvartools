# Plotting Forecasts of a List of Bayesian VAR Models

A plot function for objects of class 'expandwindbvarprdlist'.

## Usage

``` r
# S3 method for class 'expandwindbvarprdlist'
plot(x, n_pre = NULL, ci = 0.95, ...)
```

## Arguments

- x:

  an object of class 'expandwindbvarprdlist', usually, a result of a
  call to
  [`predict.expandingwindow`](https://franzmohr.github.io/bvartools/reference/predict.expandingwindow.md).

- n_pre:

  number of plotted observations that precede the forecasts. If `NULL`
  (default), all available observations will be plotted.

- ci:

  interval used to calculate the credible bands of the forecasts.

- ...:

  further graphical parameters, which are passed on to
  [`plot.ts`](https://rdrr.io/r/stats/plot.ts.html).
