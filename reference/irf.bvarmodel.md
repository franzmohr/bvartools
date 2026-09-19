# Impulse Response Function

Computes the impulse response coefficients for an object of class
'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
irf(
  x,
  impulse = NULL,
  response = NULL,
  n_ahead = 5,
  ci = 0.95,
  shock = 1,
  type = "feir",
  cumulative = FALSE,
  keep_draws = FALSE,
  period = NULL,
  impact = NULL,
  ...
)
```

## Arguments

- x:

  an object of class 'bvarmodel'.

- impulse:

  name of the impulse variable.

- response:

  name of the response variable.

- n_ahead:

  number of steps ahead. Zero is allowed and returns the impact response
  alone.

- ci:

  a numeric between 0 and 1 specifying the probability mass covered by
  the credible intervals. Defaults to 0.95.

- shock:

  size of the shock. For `"oir"`, `"gir"` and `"sgir"` it is counted in
  standard deviations of the shock, so the default of 1 is a shock of
  one standard deviation. For `"feir"` and `"sir"` it is counted in
  units of the reduced form or structural error of the impulse variable,
  and for `"sign"` and `"custom"` it rescales the columns of the impact
  matrix. `"sd"` and `"nsd"` are a positive and a negative shock of one
  standard deviation and are not available for `"sign"` and `"custom"`.

- type:

  type of the impulse response. Possible choices are forecast error
  `"feir"` (default), orthogonalised `"oir"`, structural `"sir"`,
  generalised `"gir"`, structural generalised `"sgir"`, sign restricted
  `"sign"` and `"custom"` impulse responses. For a structural model only
  `"sir"` and `"sgir"` are available; the other five require a
  non-structural model. See 'Details'.

- cumulative:

  logical specifying whether a cumulative IRF should be calculated.

- keep_draws:

  logical specifying whether the function should return all draws of the
  posterior impulse response function. Defaults to `FALSE` so that the
  median and the credible intervals of the posterior draws are returned.

- period:

  integer. Index of the period, for which the IR should be generated.
  Only used for TVP or SV models. Default is `NULL`, so that the
  posterior draws of the last time period are used.

- impact:

  the impact matrix of a `"custom"` impulse response, either a single
  \\K \times K\\ matrix that identifies every posterior draw the same
  way, or a list of such matrices with one entry per draw. Ignored for
  every other value of `type`.

- ...:

  further arguments passed to or from other methods.

## Value

A time-series object of class 'bvarirf' running from period 0 to
`n_ahead`, with the lower bound, the median and the upper bound of the
credible band of the response in three columns named after their
quantiles, e.g. `"2.5%"`, `"50%"` and `"97.5%"` for `ci = .95`. If
`keep_draws = TRUE`, a matrix of class 'bvarirf' with one row per draw
and one column per period instead.

## Details

The function produces different types of impulse responses for the VAR
model \$\$A_0 y_t = \sum\_{i = 1}^{p} A\_{i} y\_{t-i} + u_t,\$\$ with
\\u_t \sim N(0, \Sigma)\\.

Forecast error impulse responses \\\Phi_i\\ are obtained by recursions
\$\$\Phi_i = \sum\_{j = 1}^{i} \Phi\_{i-j} A_j, i = 1, 2,...,h\$\$ with
\\\Phi_0 = I_K\\.

Orthogonalised impulse responses \\\Theta^o_i\\ are calculated as
\\\Theta^o_i = \Phi_i P\\, where P is the lower triangular Choleski
decomposition of \\\Sigma\\, so that they are the responses to
orthogonalised shocks of one standard deviation.

Structural impulse responses \\\Theta^s_i\\ are calculated as
\\\Theta^s_i = \Phi_i A_0^{-1}\\.

For a structural model the posterior draws describe the structural form:
the coefficients are the \\A_i\\ of the equation above and \\\Sigma\\ is
the covariance of the structural errors. The recursion for \\\Phi_i\\
needs the reduced form, i.e. \\A_0^{-1} A_i\\ and \\A_0^{-1} \Sigma
A_0^{-1\prime}\\, which only `"sir"` and `"sgir"` form. The other types
are therefore not available for a structural model, and the two
structural types not for any other.

Sign restricted impulse responses \\\Theta^r_i\\ are calculated as
\\\Theta^r_i = \Phi_i P Q\\, where \\P\\ is the lower triangular
Choleski decomposition of \\\Sigma\\ and \\Q\\ the rotation that
[`add_sign_restrictions`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md)
accepted for that draw. Draws for which no admissible rotation was found
are left out, so the responses cover fewer draws than the posterior
holds and the credible interval is one over the set of models the
restrictions admit rather than over a single identified model.

Custom impulse responses \\\Theta^c_i\\ are calculated as \\\Theta^c_i =
\Phi_i P\\, where \\P\\ is the matrix supplied in argument `impact`.
This is the route by which an identification that the package does not
derive itself reaches the recursion; `shock` rescales the result but
nothing normalises the columns of \\P\\, so an impact matrix that means
to deliver unit shocks has to arrive that way.

(Structural) Generalised impulse responses to a shock to variable \\j\\
are calculated as \\\Theta^g\_{i} = \sigma\_{jj}^{-1/2} \Phi_i A_0^{-1}
\Sigma e_j\\, where \\\sigma\_{jj}\\ is the \\j\\th diagonal element of
\\\Sigma\\, the variance of the error of the impulse variable, and
\\e_j\\ is a selection vector containing one in its \\j\\th element and
zero otherwise. They are therefore the responses to a shock of one
standard deviation (Pesaran and Shin, 1998). If the `"bvarmodel"` object
does not contain draws of \\A_0\\, it is assumed to be an identity
matrix.

## References

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis
in linear multivariate models. *Economics Letters, 58*, 17-29.

## See also

Other post-estimation analysis:
[`add_sign_restrictions.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md),
[`fevd.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md),
[`fevd.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md),
[`irf.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md),
[`predict.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md),
[`spillover.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md),
[`spillover.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md),
[`vec_to_var.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Obtain IR
ir <- irf(model, impulse = "invest", response = "cons")

```
