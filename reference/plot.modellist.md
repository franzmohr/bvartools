# Plotting Draws of a Bayesian Time Series Models

A plot function for objects of class 'modellist'.

## Usage

``` r
# S3 method for class 'modellist'
plot(x, ...)
```

## Arguments

- x:

  an object of class 'modellist'.

- ...:

  arguments passed forward to other methods.

## Value

`x`, invisibly. The function is called for its side effect, one plot per
model in the list, drawn by the plot method of that model.
