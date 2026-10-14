# Check a volatility break and place it in the sample

Validates the volatility break of Lenza and Primiceri (2022) and returns
the two fields a model carries for it, or `NULL` when no break was asked
for.

## Usage

``` r
check_volatility_break(volatility_break, volatility_scales, y, algorithm)
```

## Arguments

- volatility_break:

  the first period of the break as a date of the data, `c(year, period)`
  such as `c(2020, 3)`, or a single time such as `2020`. `NULL` for no
  break.

- volatility_scales:

  a whole number of at least one, how many periods from the break on get
  a scale of their own before the excess decays.

- y:

  the endogenous series of the estimation sample, whose
  [`time`](https://rdrr.io/r/stats/time.html) the date is looked up in.

- algorithm:

  the sampler the model would be estimated with.

## Value

A named list with the integer elements `volatility_break` and
`volatility_scales`, or `NULL`.

## Details

[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
and
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
call this, and it is exported so that a package building its own models
– the sub-models of a global model in bgvars among them – validates a
break the same way rather than repeating the arithmetic. The break is
given as a date and kept as the period of `y` it falls on, counted from
one, which is what `/model/volatility_break` of a model file holds and
what the sampler reads.

**It is refused for every algorithm but the two constant coefficient
Wishart samplers**, `VarNormalWishart` and `VecNormalWishart`, which are
the ones that scale the error covariance from a known period on. A break
outside the estimation sample is refused as well: the periods after it
are what estimate its size.

## References

Lenza, M. and Primiceri, G. E. (2022). How to estimate a vector
autoregression after March 2020. *Journal of Applied Econometrics*,
37(4), 688–699.

## See also

[`volatility_break_prior`](https://franzmohr.github.io/bvartools/reference/volatility_break_prior.md)
