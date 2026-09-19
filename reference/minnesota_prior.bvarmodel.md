# Minnesota Prior

Calculates the Minnesota prior for a VAR model.

## Usage

``` r
# S3 method for class 'bvarmodel'
minnesota_prior(
  object,
  kappa1 = 2,
  kappa2 = 0.5,
  kappa3 = NULL,
  kappa4 = 5,
  max_var = NULL,
  coint_var = FALSE,
  sigma = "AR",
  ...
)
```

## Arguments

- object:

  an object of class 'bvarmodel', usually, a result of a call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md).

- kappa1:

  a numeric specifying the prior variance of coefficients that
  correspond to own lags of endogenous variables. See 'Details'.

- kappa2:

  a numeric specifying the size of the prior variance of endogenous
  variables, which do not correspond to own lags. See 'Details'.

- kappa3:

  a numeric specifying the size of the prior variance of
  non-deterministic exogenous variables. Default is `NULL`, which
  indicates that the formula for the calculation of the prior variance
  of deterministic terms is used for all exogenous variables. See
  'Details'.

- kappa4:

  a numeric specifying the size of the prior variance of deterministic
  terms relative to argument `kappa1`.

- max_var:

  a positive numeric specifying the maximum prior variance that is
  allowed for coefficients of non-deterministic variables. If `NULL`
  (default), the prior variances are not limited.

- coint_var:

  a logical specifying whether the model is a cointegrated VAR model,
  for which the prior means of first own lags should be set to one.

- sigma:

  either `"AR"` (default) or `"VAR"` indicating that the variances of
  the endogenous variables \\\sigma^2\\ are calculated based on a
  univariate AR regression or a least squares estimate of the VAR form,
  respectively. In both cases all deterministic variables are used in
  the regressions, if they appear in the model.

- ...:

  further arguments passed to or from other methods.

## Value

A list containing a matrix of prior means and the precision matrix of
the coefficients and the inverse variance-covariance matrix of the error
term, which was obtained by an LS estimation.

## Details

The function calculates the Minnesota prior of a VAR model. For the
endogenous variable \\i\\ the prior variance of the \\l\\th lag of
regressor \\j\\ is obtained as \$\$ \frac{\kappa\_{1}}{l^2} \textrm{ for
own lags of endogenous variables,}\$\$ \$\$ \frac{\kappa\_{1}
\kappa\_{2}}{l^2} \frac{\sigma\_{i}^2}{\sigma\_{j}^2} \textrm{ for
endogenous variables other than own lags,}\$\$ \$\$ \frac{\kappa\_{1}
\kappa\_{3}}{(l+1)^2} \frac{\sigma\_{i}^2}{\sigma\_{j}^2} \textrm{ for
unmodelled exogenous variables,}\$\$ \$\$ \kappa\_{1} \kappa\_{4}
\sigma\_{i}^2 \textrm{ for deterministic terms,}\$\$ where
\\\sigma\_{i}\\ is the residual standard deviation of variable \\i\\ of
an unrestricted LS estimate. For exogenous variables \\\sigma\_{i}\\ is
the sample standard deviation. In case structural parameters are
estimated, the formula \\\kappa\_{1} \kappa\_{2}
\frac{\sigma\_{i}^2}{\sigma\_{j}^2}\\ is used. If \\kappa\_{3}\\ is not
provided, prior variances are calculated in the same way as for
deterministic terms. If the model does not contain exogenous variables,
argument `kappa3` will be ignored.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). *Bayesian
Econometric Methods* (2nd ed.). Cambridge: University Press.

Lütkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model input
object <- create_bvarmodel(e6, p = 1,
                           deterministic = "both",
                           seasonal = TRUE)

# Obtain Minnesota prior
prior <- minnesota_prior(object)
```
