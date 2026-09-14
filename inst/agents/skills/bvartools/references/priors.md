# Priors: what `add_priors()` needs

Every `r` block on this page runs in the package's test suite, top to bottom in
one session, so the shapes and errors asserted with `stopifnot()` are what the
installed version does. Draw counts are kept small; the priors do not depend on
them.

`add_priors()` takes named lists and **has no default hyperparameters**. Every
value a model needs must be given, in the list it belongs to:

| Argument | For | Holds |
| --- | --- | --- |
| `coef` | every model | the prior of the coefficients |
| `sigma` | every model | the prior of the error term, by `error` |
| `coint` | VEC models | the prior of the cointegration space |
| `varsel` | `varsel = "ssvs"` or `"bvs"` | the variable selection prior; an error otherwise |

A missing required element stops with a message naming it, and so does an
element that is **not recognised**, in any of the four lists. Variances are given
as **precisions**: `v_i = 0` is uninformative and a larger `v_i` shrinks harder.
The exception is `coef` of a model with `tvp = TRUE`, where `v_i` and `v_i_det`
are the prior precision of the coefficients before the sample and must be
positive: the sampler integrates that state out of the first period's prior,
which takes the inverse of the precision, and refuses a zero one.

The page to read for the full list is the method's, not the generic's:
`?add_priors.bvarmodel` and `?add_priors.bvecmodel` differ.

## `coef`

| Element | Meaning |
| --- | --- |
| `v_i` | Prior precision of the coefficients. Required unless `minnesota` is given, and required with it for `error = "gamma+covar"` or `"sv+covar"`, where it is also the precision of the covariance coefficients |
| `v_i_det` | Prior precision of the deterministic terms. Falls back to `v_i` |
| `const` | Prior mean of the intercepts: a number, one per variable, `"mean"` or `"first"` |
| `coint_var` | VAR only. `TRUE` sets the prior mean of the first own lag to one |
| `minnesota` | A list with `kappa1`, `kappa2`, `kappa4`, and `kappa3` if there are exogenous variables. Replaces `v_i` and `v_i_det` for the coefficients |
| `max_var` | Cap on the Minnesota prior variances of the non-deterministic coefficients |
| `shape`, `rate` | Required for `tvp = TRUE`: the gamma prior of the state variances |
| `rate_det` | Rate for the state variances of deterministic terms. Falls back to `rate` |

The coefficient vector holds all lag coefficients first, then the exogenous
ones, then the deterministic ones, so `v_i_det` lands at the end:

```r
library(bvartools)
set.seed(1)

data("e1")
e1 <- diff(log(e1)) * 100
k <- ncol(e1)

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 50)
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 0.1),
                    sigma = list(df = "k + 2", scale = 1))

precision <- diag(model$priors$a$v_inv)
n_lags <- k * k * 2                      # two lags of three variables, per equation
stopifnot(length(precision) == n_lags + k,
          all(precision[seq_len(n_lags)] == 1),
          all(precision[n_lags + seq_len(k)] == 0.1),
          model$priors$u_sigma$df == k + 2)
```

A misspelt element is an error, not a silently ignored setting:

```r
misspelt <- tryCatch(add_priors(model, coef = list(v_i = 1),
                                sigma = list(df = "k", scal = 1)),
                     error = function(e) conditionMessage(e))
stopifnot(grepl("not recognised", misspelt))
```

A Minnesota prior shrinks more distant lags harder:

```r
minn <- add_priors(model,
                   coef = list(minnesota = list(kappa1 = 0.2, kappa2 = 0.5, kappa4 = 100)),
                   sigma = list(df = "k", scale = 1))

variance <- 1 / diag(minn$priors$a$v_inv)
first_lag <- variance[seq_len(k * k)]
second_lag <- variance[k * k + seq_len(k * k)]
stopifnot(all(second_lag < first_lag))
```

## `sigma`, by `error`

`error` is fixed in `create_bvarmodel()` or `create_bvecmodel()`, and decides what
`sigma` must contain:

| `error` | `sigma` must contain | Prior |
| --- | --- | --- |
| `"wishart"` (default) | `df`, `scale` | Inverse Wishart on the covariance matrix |
| `"gamma"`, `"gamma+covar"` | `shape`, `rate` | Gamma on the error precisions |
| `"sv"`, `"sv+covar"` | `mu`, `v_i`, `shape`, `rate`, `state_variance`, `offset` | Stochastic volatility |
| `"ald"` (VAR only) | `shape`, `rate` | Scales of the asymmetric Laplace distribution |

- `df` and a gamma `shape` accept a number or an expression in `k`, the number
  of endogenous variables, such as `"k"` or `"k + 3"`, and are stored as given,
  for a VEC as for a VAR. The Wishart VEC samplers add the rank `r` to the
  posterior degrees of freedom themselves.
- `+covar` estimates the error covariances through a triangular decomposition;
  without it the error term is diagonal.
- For stochastic volatility, `shape` and `rate` are the prior of the variance of
  the log-volatility state equation, `mu` and `v_i` the prior of its initial
  state, `state_variance` the starting value of that variance, and `offset` the
  constant added before taking logs of squared errors.

A structural model cannot use the Wishart prior: `A_0` is not identified against
an unrestricted covariance. Use a diagonal error term:

```r
refused <- tryCatch(create_bvarmodel(e1, p = 1, structural = TRUE, error = "wishart",
                                     iterations = 100, burnin = 50),
                    error = function(e) "refused")
stopifnot(identical(refused, "refused"))

structural <- create_bvarmodel(e1, p = 1, deterministic = "const",
                               structural = TRUE, error = "gamma",
                               iterations = 100, burnin = 50)
structural <- add_priors(structural,
                         coef = list(v_i = 1),
                         sigma = list(shape = 3, rate = 1))
stopifnot(identical(structural$priors$u_sigma$type, "gamma"))
```

## `varsel`

| Element | For | Meaning |
| --- | --- | --- |
| `inprior` | both | Prior inclusion probability, between 0 and 1. Required |
| `tau` | SSVS | Prior standard deviations of excluded and included coefficients, `c(tau0, tau1)` |
| `semiautomatic` | SSVS | Two factors scaling least squares standard errors into `tau0` and `tau1`, the approach of George et al. (2008) |
| `exclude_det` | both | `TRUE` keeps deterministic terms out of the selection |
| `covar` | both | `TRUE` also selects over the error covariances |
| `minnesota` | both | Four numbers for Minnesota-like inclusion probabilities |

SSVS needs `tau` or `semiautomatic`; BVS needs only `inprior`:

```r
bvs <- create_bvarmodel(e1, p = 2, deterministic = "const", varsel = "bvs",
                        iterations = 100, burnin = 50)
bvs <- add_priors(bvs,
                  coef = list(v_i = 1, v_i_det = 0.1),
                  sigma = list(df = "k", scale = 1),
                  varsel = list(inprior = 0.5, exclude_det = TRUE))
stopifnot(all(c("inprior", "include") %in% names(bvs$priors$a)))

ssvs <- create_bvarmodel(e1, p = 2, deterministic = "const", varsel = "ssvs",
                         iterations = 100, burnin = 50)
ssvs <- add_priors(ssvs,
                   coef = list(v_i = 0),
                   sigma = list(df = "k", scale = 1),
                   varsel = list(inprior = 0.5, semiautomatic = c(0.01, 10)))
stopifnot(all(c("tau0", "tau1") %in% names(ssvs$priors$a)))
```

Refused, each with a message saying why:

- SSVS with stochastic volatility, and SSVS with time-varying parameters.
- SSVS for a quantile VAR (`error = "ald"`); use BVS there.
- Selection over the covariances (`covar = TRUE`) with a Wishart prior.
- A `varsel` list for a model created with `varsel = "none"`.

## `coint`, for VEC models

| Element | For | Meaning |
| --- | --- | --- |
| `v_i` | constant cointegration | Shrinkage of the cointegration space prior, or `"ml"`. `0` gives a uniform prior on the space |
| `p_tau_i` | constant cointegration | Inverse of the central location matrix of the space: diagonal elements, a full matrix, or `"ml"` to centre it on Johansen's estimate |
| `weight` | `p_tau_i = "ml"` | Weight of the estimate, in units of the information in the sample. Defaults to 1 |
| `rho` | time-varying cointegration | Autocorrelation of the state equation, below one. Required |
| `rho_min`, `rho_max` | time-varying cointegration | Support of a uniform prior on `rho`, both or neither |

```r
data("e6")
e6 <- e6 * 100

vec <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
                        iterations = 100, burnin = 50)
vec <- add_priors(vec,
                  coef = list(v_i = 0, v_i_det = 0),
                  coint = list(v_i = 0, p_tau_i = 1),
                  sigma = list(df = "k", scale = 1))

stopifnot(identical(vec$priors$beta$type, "cointspace"),
          vec$priors$u_sigma$df == ncol(e6))   # "k", as given
```

A prior centred on the maximum likelihood estimate goes with starting values
from it, as in `vignette("bvec", package = "bvartools")`:

```r
vec_ml <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
                           iterations = 100, burnin = 50)
vec_ml <- add_priors(vec_ml,
                     coef = list(v_i = 0, v_i_det = 0),
                     coint = list(v_i = "ml", p_tau_i = "ml", weight = 1),
                     sigma = list(df = "k", scale = 1))
vec_ml <- add_initial_values(vec_ml, method = "maxlik")
stopifnot(identical(vec_ml$priors$beta$type, "cointspace"))
```

## Checking what was set

After `add_priors()`, look at the result rather than trusting the call:

```r
str(model$priors, max.level = 2)
```
