# Thinning Posterior Draws

Thins the MCMC posterior draws of the elements in an object of class
'modellist'.

## Usage

``` r
# S3 method for class 'modellist'
thin(x, thin = 10, ...)
```

## Arguments

- x:

  an object of class 'modellist'.

- thin:

  an integer specifying the thinning interval between successive values
  of posterior draws.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'modellist'.
