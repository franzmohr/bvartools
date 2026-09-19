# Austrian sub-model of a global VAR

The data set contains quarterly time series for Austria, the
corresponding foreign (star) variables and global commodity prices from
1979Q2 to 2023Q3. It was produced from data set `gvar2023` of package
`bgvars`, which contains the GVAR database of Mohaddes and Raissi
(2024), as the Austrian sub-model of a global VAR model of all 33
countries of the database.

## Usage

``` r
data("at_macrodata")
```

## Format

A named list with two elements, which are quarterly time-series objects
from 1979Q2 to 2023Q3 with 178 rows each:

- `domestic`:

  a time-series object with 6 columns, the domestic (endogenous)
  variables:

  `y`

  :   log real GDP.

  `Dp`

  :   rate of inflation, the quarterly change in the log CPI.

  `eq`

  :   log real equity prices.

  `ep`

  :   log exchange rate against the US dollar, deflated by the CPI.

  `r`

  :   short-term interest rate, \\0.25 \ln(1 + R^S / 100)\\, where
      \\R^S\\ is the annual short-term rate in percent.

  `lr`

  :   long-term interest rate, \\0.25 \ln(1 + R^L / 100)\\, where
      \\R^L\\ is the annual long-term rate in percent.

- `foreign`:

  a time-series object with 9 columns, the foreign and global (weakly
  exogenous) variables:

  `y.s`, `Dp.s`, `eq.s`, `ep.s`, `r.s`, `lr.s`

  :   foreign counterparts of the domestic variables.

  `poil`

  :   log oil prices.

  `pmat`

  :   log agricultural raw material prices.

  `pmetal`

  :   log metals prices.

## Details

The foreign variables are trade weighted averages of the series of the
other countries of the database for which the respective variable is
available. The weights are constant over the whole sample and equal the
shares of the countries in Austria's trade summed over the years 2014 to
2016.

In a sub-model of a global VAR the domestic variables are endogenous,
and the foreign and global variables are weakly exogenous. The two
elements of the list can therefore be passed to arguments `data` and
`exogen` of
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
or
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
e.g.
`create_bvarmodel(data = at_macrodata$domestic, exogen = at_macrodata$foreign)`.

## References

Dees, S., di Mauro, F., Pesaran, M. H., & Smith, L. V. (2007). Exploring
the international linkages of the euro area: A global VAR analysis.
*Journal of Applied Econometrics, 22*(1), 1–38.
[doi:10.1002/jae.932](https://doi.org/10.1002/jae.932)

Mohaddes, K., & Raissi, M. (2024). *Compilation, revision and updating
of the global VAR (GVAR) database, 1979Q2–2023Q3* (mimeo). University of
Cambridge: Judge Business School.
