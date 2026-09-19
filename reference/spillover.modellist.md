# Spillover Index

Produces the connectedness measures of Diebold and Yilmaz (2012) for
every model in an object of class 'modellist'.

## Usage

``` r
# S3 method for class 'modellist'
spillover(object, ...)
```

## Arguments

- object:

  an object of class 'modellist'.

- ...:

  arguments passed forward to
  [`spillover.bvarmodel`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md).

## Value

A list of objects of class 'bvarspillover'. Elements belonging to models
whose estimation failed are `NULL`.
