# The Discounted Models

Whether a model is one of the two discounted models, and what a
discounted model cannot be. Both are exported for the packages that
build a model of their own on this one and have to make the same
distinction.

## Usage

``` r
is_discount_model(x)

check_discount_specification(
  k,
  error = "wishart",
  varsel = "none",
  structural = FALSE,
  burnin = 0,
  thin = 1,
  delta_beta = 1,
  delta_sigma = 1
)
```

## Arguments

- x:

  a model, a specification list – the `model` element of one – or the
  name of an algorithm.

- k:

  the number of endogenous variables.

- error, varsel, structural, burnin, thin:

  the corresponding arguments of
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md).

- delta_beta, delta_sigma:

  the two discount factors.

## Value

`is_discount_model` returns a single logical.
`check_discount_specification` returns `NULL` invisibly.

## Details

`VarTvpDiscount` and `VecTvpDiscount` are the matrix normal dynamic
linear model of West & Harrison (1997, ch. 16) with the discounted
Wishart of Uhlig (1997). Their posterior is closed form – one pass over
the sample, no chain, no random numbers consumed – which is what every
refusal below follows from. See
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
for what they are and how they are set up.

`check_discount_specification` raises an error for a specification a
discounted model cannot carry and returns nothing otherwise. Each
refusal is one BayesTS would raise as well, against a file that had
already been written; raising it here means one error for a grid of
models rather than one per file.

## References

Uhlig, H. (1997). Bayesian vector autoregressions with stochastic
volatility. *Econometrica, 65*(1), 59–73.
[doi:10.2307/2171813](https://doi.org/10.2307/2171813)

West, M., & Harrison, J. (1997). *Bayesian forecasting and dynamic
models* (2nd ed.). New York: Springer.

## Examples

``` r

data("e6")
model <- create_bvecmodel(e6 * 100, p = 2, r = 1, const = "unrestricted",
                          algorithm = "discount", delta_beta = 0.98,
                          iterations = 10, burnin = 0, thin = 1)
is_discount_model(model)
#> [1] TRUE
is_discount_model("VecNormalWishart")
#> [1] FALSE
```
