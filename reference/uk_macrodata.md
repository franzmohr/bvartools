# UK interest and inflation rate data

The data set contains quarterly time series for two UK interest rates
and inflation from 1957Q4 to 2009Q3. It was produced from the
supplementary material to Koop, León-González and Strachan (2011).

## Usage

``` r
data("uk_macrodata")
```

## Format

A named time-series object with 210 rows and 3 variables:

- rs:

  short term rate: Treasury bills. This should be series 'FITB_PA' of
  the IMF's International Financial Statistics database.

- rs:

  long term rate: Government bonds. This should be series 'FIGB_PA' of
  the IMF's International Financial Statistics database.

- Dp:

  inflation: Quarterly change in the log CPI annualized by multiplying
  the change by 400.

## References

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference
in a time varying cointegration model. *Journal of Econometrics,
165*(2), 210–220.
[doi:10.1016/j.jeconom.2011.07.007](https://doi.org/10.1016/j.jeconom.2011.07.007)
