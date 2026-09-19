# US macroeconomic data

The data set contains quarterly time series for the US CPI inflation
rate, unemployment rate, and Fed Funds rate from 1959Q2 to 2007Q4. It
was produced from file "US_macrodata.csv" of the data sets associated
with Chan, Koop, Poirier and Tobias (2019). Raw data are available at
<https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_1/US_macrodata.csv>.

## Usage

``` r
data("us_macrodata")
```

## Format

A named time-series object with 195 rows and 3 variables:

- Dp:

  CPI inflation rate.

- u:

  unemployment rate.

- r:

  Fed Funds rate.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.
