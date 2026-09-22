# Apply Transformations

Transforms the columns of a time-series object using the seven
transformation codes of FRED-MD/FRED-QD.

## Usage

``` r
transform_variables(x, code)
```

## Arguments

- x:

  a time-series object.

- code:

  a named integer vector of transformation codes, one of `1:7` – see
  'Details'. Names must match columns of `x`; a column with no entry is
  left untransformed (code `1`). For a single series `code` may also be
  one unnamed code.

## Value

A time-series object of the same shape as `x`: a vector series for a
vector series, and a matrix of series with the same columns otherwise.

## Details

The seven codes, applied column by column, are

- 1:

  no transformation, \\x_t\\

- 2:

  first difference, \\\Delta x_t\\

- 3:

  second difference, \\\Delta^2 x_t\\

- 4:

  logarithm, \\\log(x_t)\\

- 5:

  first difference of the logarithm, \\\Delta \log(x_t)\\

- 6:

  second difference of the logarithm, \\\Delta^2 \log(x_t)\\

- 7:

  first difference of the growth rate, \\\Delta(x_t / x\_{t - 1} - 1)\\

Codes `4:6` need a strictly positive series.

A transformation loses as many leading observations as it reaches back –
one for codes `2` and `5`, two for `3`, `6` and `7`, the last being the
difference of a growth rate that is itself taken over one period –
rather than shortening the series: those leading periods become `NA`, so
every column keeps the time index of `x` and can still be combined with
[`ts.intersect`](https://rdrr.io/r/stats/ts.union.html) or passed to
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

## References

McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for
macroeconomic research. *Journal of Business & Economic Statistics,
34*(4), 574–589.

McCracken, M. W., & Ng, S. (2021). FRED-QD: A quarterly database for
macroeconomic research. *Federal Reserve Bank of St. Louis Review,
103*(1), 1–44.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

data("us_macrodata")

code <- c(Dp = 2, u = 1, r = 2)
x <- transform_variables(us_macrodata, code)
```
