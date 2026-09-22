# Test for Time Variation in a VAR or VEC Model

Computes the Bayes factors of Chan (2018) for time variation in each
coefficient, each covariance coefficient and, under stochastic
volatility, each log-volatility of a VAR or VEC model with time varying
parameters, from the draws of a single estimation of that model.

## Usage

``` r
# S3 method for class 'bvarmodel'
time_variation_test(object, joint = TRUE, batches = 20, ...)

# S3 method for class 'bvecmodel'
time_variation_test(object, joint = TRUE, batches = 20, ...)
```

## Arguments

- object:

  an object of class `"bvarmodel"` with `tvp = TRUE` and `error = "sv"`,
  `"sv+covar"`, `"gamma"` or `"gamma+covar"`, or of class `"bvecmodel"`
  with the same, whose priors were set with `coef$omega_v` or, under
  stochastic volatility, `sigma$omega_v` in
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
  and whose posterior was drawn by
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- joint:

  logical. Should the Bayes factor for the time variation of every state
  of a block at once be reported as well? Default is `TRUE`.

- batches:

  integer. The number of batches the draws are split into for the
  numerical standard errors. Default is 20.

- ...:

  further arguments, which are ignored.

## Value

A data frame of class `"bvartimevar"` with one row per state and, if
`joint = TRUE`, one per block, and the columns

- `block`:

  `"coefficients"`, `"covariances"` or `"volatilities"`.

- `equation`:

  the endogenous variable of the equation the state belongs to; `NA` for
  the joint rows.

- `term`:

  the regressor of a coefficient, the variable whose error a covariance
  coefficient loads on, `"log-volatility"`, or `"(joint)"`.

- `log_bf`:

  the log Bayes factor in favour of time variation.

- `nse`:

  its numerical standard error.

## Details

Under the prior `omega_v` a random walk \\x_t = x\_{t-1} + v_t\\, \\v_t
\sim N(0, \omega^2)\\, is estimated in the non-centred form \\x_t =
x_0 + \omega \tilde{x}\_t\\ of Frühwirth-Schnatter and Wagner (2010),
with \\\tilde{x}\_t\\ a standard random walk and the prior \\\omega \sim
N(0, V\_\omega)\\. The state does not move exactly when \\\omega = 0\\,
which is a point inside the prior, so the Bayes factor of the model in
which the state varies against the one in which it is constant is the
Savage-Dickey density ratio \$\$BF = \frac{p(\omega = 0)}{p(\omega = 0
\| y)}.\$\$ The numerator is the density of the prior at zero. The
denominator is estimated by the average over the draws of the density at
zero of the conditional posterior of \\\omega\\, which the sampler
stores for every draw in `omega_log_zero` (Chan 2018). A positive log
Bayes factor favours time variation. On the scale of Kass and Raftery
(1995), values of the log Bayes factor between 1 and 3 are positive
evidence, between 3 and 5 strong evidence and above 5 very strong
evidence, and the same values with a negative sign are evidence for a
constant state.

The joint Bayes factor of a block compares the model in which every
state of the block varies with the model in which none does. It is not a
test of whether at least one state varies: every state that is constant
costs the joint model about the same as one of its own Bayes factors
against time variation, so a block in which one state moves and many do
not can come out against time variation as a whole. Conversely, the
joint Bayes factor can be far larger than the individual ones: each
state may have posterior mass near zero on its own while little of the
joint posterior sits where all of them are near zero at once. Chan
(2018, Table 2 and footnote 8) reports log Bayes factors of about 3 for
each of two volatilities of Italian inflation and of 235 for the two
together.

The estimate is least precise where the Bayes factor is large, since the
ordinates are then an average of small numbers (Chan 2018, section 2.3).
The numerical standard error of each log Bayes factor is obtained by the
delta method from batch means, and says how much of the value is owed to
the length of the chain. Where it is large, more draws are needed before
the value is read as more than its sign.

Each block chooses its prior separately in
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md),
and only the blocks estimated under `omega_v` are reported.

For a VEC model the coefficients are the loadings, labelled by the error
correction term they load on (`ect1`, `ect2`, ...), followed by the
other coefficients. The cointegration space itself follows a state
equation with a fixed variance, which has no prior to test against and
is not reported.

## References

Chan, J. C. C. (2018). Specification tests for time-varying parameter
models with stochastic volatility. *Econometric Reviews, 37*(8),
807–823.
[doi:10.1080/07474938.2016.1167948](https://doi.org/10.1080/07474938.2016.1167948)

Frühwirth-Schnatter, S., & Wagner, H. (2010). Stochastic model
specification search for Gaussian and partial non-Gaussian state space
models. *Journal of Econometrics, 154*(1), 85–100.
[doi:10.1016/j.jeconom.2009.07.003](https://doi.org/10.1016/j.jeconom.2009.07.003)

Kass, R. E., & Raftery, A. E. (1995). Bayes factors. *Journal of the
American Statistical Association, 90*(430), 773–795.
[doi:10.1080/01621459.1995.10476572](https://doi.org/10.1080/01621459.1995.10476572)

## Examples

``` r
# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model
model <- create_bvarmodel(e1, p = 1, deterministic = "const", tvp = TRUE,
                          error = "sv", iterations = 100, burnin = 50)

# The non-centred prior on the coefficients and the log-volatilities
model <- add_priors(model,
                    coef = list(v_i = 1 / 10, omega_v = 0.001),
                    sigma = list(mu = 0, v_i = 1 / 100, omega_v = 0.1,
                                 state_variance = 0.05, offset = 1e-4))
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# Bayes factors for time variation
time_variation_test(model)
#> Bayes factors for time variation (Savage-Dickey, Chan 2018)
#> log BF > 0 favours time variation; |log BF| > 3 is strong evidence
#> 
#> Coefficients:
#>  equation term      log BF NSE 
#>  invest   invest.l1 -0.09  0.03
#>  income   invest.l1 -1.04  0.10
#>  cons     invest.l1 -0.84  0.08
#>  invest   income.l1 -0.05  0.01
#>  income   income.l1 -0.17  0.08
#>  cons     income.l1 0.33   0.11
#>  invest   cons.l1   -0.03  0.02
#>  income   cons.l1   0.70   0.17
#>  cons     cons.l1   0.15   0.08
#>  invest   const     -0.01  0.00
#>  income   const     0.34   0.08
#>  cons     const     0.16   0.05
#>           (joint)   0.71   0.37
#> 
#> Volatilities:
#>  equation term           log BF NSE 
#>  invest   log-volatility 6.69   0.79
#>  income   log-volatility 3.87   0.80
#>  cons     log-volatility -0.12  0.51
#>           (joint)        24.87  0.98
#> 
```
