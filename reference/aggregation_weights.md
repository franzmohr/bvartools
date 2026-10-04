# Weights of a Temporal Aggregate

The weights with which a series observed at a lower frequency aggregates
the high-frequency periods of a model, for argument `aggregate` of
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

## Usage

``` r
aggregation_weights(type = c("average", "sum", "growth"), n = 3)
```

## Arguments

- type:

  a character naming the aggregation: `"average"` for a series that is
  the mean of the periods it covers, such as a quarterly unemployment
  rate in a monthly model, `"sum"` for a flow that is their total, and
  `"growth"` for the growth rate of an average, following Mariano and
  Murasawa (2003). See 'Details'.

- n:

  an integer, the number of high-frequency periods one low-frequency
  period covers: 3 for quarters of months, 12 for years of months, 4 for
  years of quarters.

## Value

A numeric vector of weights, oldest period first.

## Details

An observation of the low-frequency series in period \\t\\ is \$\$x_t =
\sum\_{j = 1}^{L} w_j y\_{t - L + j},\$\$ where \\y_t\\ is the
high-frequency series the model contains and \\w_1, \dots, w_L\\ are the
weights returned, oldest first.

For `"average"` they are \\n\\ weights of \\1/n\\, and for `"sum"` \\n\\
weights of one. For `"growth"` the high-frequency series is the growth
rate – the log difference – of a series whose average the low-frequency
series is the growth rate of. Approximating the log of an average by the
average of the logs, that growth rate is a triangular combination of
\\2n - 1\\ high-frequency growth rates, which for quarters of months is
\\(1, 2, 3, 2, 1) / 3\\.

## References

Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of
business cycles based on monthly and quarterly series. *Journal of
Applied Econometrics, 18*(4), 427–443.
[doi:10.1002/jae.695](https://doi.org/10.1002/jae.695)

## See also

Other model set-up:
[`add_dummy_variables.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_dummy_variables.bvarmodel.md),
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_prior_options()`](https://franzmohr.github.io/bvartools/reference/add_prior_options.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r
aggregation_weights("growth", 3)
#> [1] 0.3333333 0.6666667 1.0000000 0.6666667 0.3333333
```
