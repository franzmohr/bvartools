# Plotting Selection Criteria

A plot function for objects of class 'selcritlist'.

## Usage

``` r
# S3 method for class 'selcritlist'
plot(x, criterion = "WAIC", ...)
```

## Arguments

- x:

  an object of class 'selcritlist', usually, a result of a call to
  [`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md).

- criterion:

  the selection criterion that should be plotted. Available choices are
  the in-sample criteria `"LL"`, `"AIC"`, `"BIC"`, `"HQ"`, `"WAIC"`
  (default), `"LOOIC"` and the out-of-sample statistics `"FE"`, `"AFE"`
  and `"RSFE"`.

- ...:

  further graphical parameters.

## Details

For in-sample criteria each model is represented by a horizontal error
bar, which is marked by the value of the criterion. The log-likelihood
has a posterior distribution, so its bar covers the credible band and is
marked by the median as well as the mean of its draws. `"WAIC"` and
`"LOOIC"` are point estimates whose bars cover a normal interval built
from their standard error. `"AIC"`, `"BIC"` and `"HQ"` are point
estimates without a standard error, so they are drawn as a single point.
The default is `"WAIC"`, for the reasons given in
[`selection_criteria`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md).
The criterion is measured on the x-axis and the models are arranged
along the y-axis, beginning with the first model of `x` at the top. This
keeps the plot readable, if `x` contains many models. The position of a
model in `x` is added at the upper end of its error bar.

The specifications of the models are used as labels of the y-axis, where
only the specifications that differ across the models are used, as in
[`print.modellist`](https://franzmohr.github.io/bvartools/reference/print.modellist.md).
The lag order is always used. The left margin is adjusted to the width
of the labels.

No title is added by default. Space above the plot is only used, if
argument `main` is specified.

For out-of-sample statistics each row of the plot corresponds to an
endogenous variable and the x-axis to the forecast horizon. Within a
forecast horizon each model is represented by an error bar, which covers
the credible band of the respective statistic and which is marked by the
median and the mean of its posterior draws. The position of a model in
`x` is added above its error bar. Since forecast errors can be positive
and negative, a horizontal reference line is added at zero for
`criterion = "FE"`.

For both types of criteria arguments `col`, `pch`, `cex`, `lwd`, `main`
and `xlab` can be used to change the appearance of the plot, where `col`
and `lwd` are recycled over the models in `x`.
