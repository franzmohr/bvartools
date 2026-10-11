# Minnesota Prior

Forwards the elements of an object of class 'modellist' forward to its
respective method.

## Usage

``` r
# S3 method for class 'expandingwindow'
minnesota_prior(object, ...)

# S3 method for class 'modellist'
minnesota_prior(object, ...)
```

## Arguments

- object:

  a list of class 'modellist'.

- ...:

  further arguments passed to or from other methods.

  The same for an object of class 'expandingwindow', each window getting
  the prior of its own sample.

## Value

A list of class 'modellist'.
