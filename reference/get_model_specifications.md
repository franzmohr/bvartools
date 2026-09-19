# Get Model Specifications

A generic function used to obtain the model specifications of a model.

## Usage

``` r
get_model_specifications(object, ...)
```

## Arguments

- object:

  an object with suitable input data passed forward to method.

- ...:

  arguments passed forward to method.

## Value

A one-row data frame with the specifications of the model. Which columns
are returned depends on the class of argument `object`:

- type:

  the type of the model.

- k:

  the number of endogenous variables.

- p:

  the lag order of the endogenous variables.

- m:

  the number of unmodelled, non-deterministic variables.

- s:

  the lag order of the unmodelled, non-deterministic variables.

- n:

  the number of deterministic terms. For objects of class 'bvecmodel'
  the columns `n_unrestricted` and `n_restricted` are returned instead.

- rank:

  the rank of the cointegration matrix. Only for objects of class
  'bvecmodel'.

- T:

  the number of observations used for estimation. Not for objects of
  class 'expandingwindow', where it differs across the estimation
  windows, and not for objects of class 'selcrit', which do not contain
  the input data.

- varsel:

  the used variable selection algorithm.

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 10, burnin = 10)

# Obtain model specifications
get_model_specifications(model)
#>   type k p m s n  T varsel
#> 1  VAR 3 2 0 0 1 89   none
```
