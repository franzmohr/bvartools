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
