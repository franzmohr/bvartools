# Transform VEC Models to VARs in Levels

The models of the estimation windows of an object of class
`'expandingwindow'` are transformed into their VAR representation in
levels, which is the form in which forecasts, forecast errors and
out-of-sample selection criteria are obtained for VEC models.

## Usage

``` r
# S3 method for class 'expandingwindow'
vec_to_var(object, ...)
```

## Arguments

- object:

  an object of class `'expandingwindow'`.

- ...:

  arguments passed forward to method.

## Value

An object of class `'expandingwindow'`.
