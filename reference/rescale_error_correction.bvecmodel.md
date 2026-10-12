# Rescale Error Correction

Puts the series in the error correction term of an object of class
'bvecmodel', and the posterior draws that belong to them, back in terms
of the input data.

## Usage

``` r
# S3 method for class 'bvecmodel'
rescale_error_correction(object, ...)
```

## Arguments

- object:

  object of class 'bvecmodel'.

- ...:

  arguments passed forward to method.

## Value

An object of class 'bvecmodel'.

## Details

The function transforms element `object$data$train$w`, the posterior
draws and the starting values back to the series of the input data,
based on the attributes `"scale"` and `"centre"` that
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
stored in `object$data$train$w`, and drops the attributes.

If \\D\\ is the diagonal matrix of scaling factors, function
[`scale_error_correction`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
replaced \\w_t\\ by \\D^{-1} w_t\\, so that the estimated error
correction term is \\\alpha \beta^{\prime} D^{-1} w_t\\. The draws of
\\\beta\\ are therefore multiplied by \\D^{-1}\\, while the draws of
\\\alpha\\ are not affected by the transformation and are carried over
unchanged.

If the series were centred on their means \\m\\, the estimated term is
\\\alpha \beta^{\prime} D^{-1} (w_t - m)\\, with \\D\\ the identity
matrix if they were not scaled. The draws of the unrestricted constant
are therefore shifted by \\-\alpha \beta^{\prime} D^{-1} m\\, draw by
draw and, for a model with time varying parameters, period by period
with the loadings and the cointegration vectors of that period. The
fitted values and the log-likelihood of every draw are unchanged by
that. The starting values of the constant are shifted in the same way.
If variable selection covered the deterministic terms, the draws of the
inclusion indicators of the constant are carried over as they are, so
they describe the constant of the centred series.

What is unchanged is the sample. Beyond it, a model with time varying
cointegration vectors cannot be carried on from the rescaled draws: the
state equation of \\\beta_t\\ was estimated for the transformed series,
and the constant absorbs \\\alpha_t \beta_t^{\prime} D^{-1} m\\ of the
last period only, so a step of \\\beta_t\\ would reintroduce the
intercept the centring removed. The function therefore sets
`object$model$ect_rescaled` for such a model, and
[`add_posterior_forecasts.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md)
and
[`add_predictive_loglik.bvecmodel`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md)
accept it only with `forecast_states = "hold"`.

## Examples

``` r

data("e6")
model <- create_bvecmodel(e6 * 100, p = 1, r = 1, const = "unrestricted",
                          iterations = 20, burnin = 10)
model <- scale_error_correction(model)
model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))
model <- add_posterior_coefficients(add_initial_values(model))

# The draws of beta on the scale of the data
model <- rescale_error_correction(model)
summary(model)
#> 
#> Bayesian VEC model with p = 1 
#> 
#> Endogenous variables: R, Dp
#> 
#> Variable: R 
#> 
#>            Mean      SD Naive SD Time-series SD     2.5%        50%   97.5% 
#> l.R    0.002064 0.00913 0.002041       0.002041 -0.01051  0.0004154 0.01586 
#> l.Dp  -0.006847 0.02949 0.006594       0.006594 -0.05517 -0.0019251 0.03592 
#> const -0.026588 0.07044 0.015751       0.004365 -0.15912 -0.0285645 0.07145 
#> 
#> Variable: Dp 
#> 
#>          Mean      SD Naive SD Time-series SD    2.5%     50%   97.5%  
#> l.R    0.3879 0.04633  0.01036        0.02386  0.3205  0.3846  0.4796 *
#> l.Dp  -1.2734 0.06936  0.01551        0.02301 -1.3893 -1.2709 -1.1501 *
#> const -1.8869 0.42137  0.09422        0.25725 -2.6194 -1.8946 -1.2703 *
#> 
#> Variance-covariance matrix:
#> 
#>          Mean      SD Naive SD Time-series SD    2.5%     50%   97.5%  
#> R_R    0.3147 0.04879  0.01091        0.01091  0.2253  0.3068  0.4059 *
#> R_Dp  -0.1489 0.07898  0.01766        0.01766 -0.2783 -0.1644 -0.0301 *
#> Dp_Dp  3.0461 0.27414  0.06130        0.06130  2.7052  2.9906  3.7230 *
#> 
```
