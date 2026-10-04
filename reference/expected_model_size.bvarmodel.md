# Expected Size of a Model

Calculates how much memory a model will take once its posterior draws
are complete, block by block, from its specification alone, so that it
can be called before a long estimation is started.

## Usage

``` r
# S3 method for class 'bvarmodel'
expected_model_size(object, chains = NULL, ...)

# S3 method for class 'bvecmodel'
expected_model_size(object, chains = NULL, ...)

# S3 method for class 'expandingwindow'
expected_model_size(object, ...)

# S3 method for class 'modellist'
expected_model_size(object, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel' or 'bvecmodel', usually a result of a
  call to
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
  or a list of such models of class 'modellist' or 'expandingwindow'.

- chains:

  the number of chains that will be simulated. If `NULL` (default), the
  value in `object$model$chains` is used, and one chain when there is
  none. See
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- ...:

  further arguments passed to or from other methods.

## Value

A data frame of class 'modelsize' with one row per block of the object
and the columns

- `model`:

  for a list of models, the name or position of the model.

- `element`:

  where the block will be in the object, such as `posterior$a$coeffs`.
  The first row of each model is everything that is there already, the
  data, the priors and the initial values.

- `step`:

  the function that adds the block, empty for the first row. The warning
  of
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  counts the first row and the blocks that function adds itself.

- `draws`, `columns`:

  the dimensions of the block.

- `bytes`:

  the size of the block in bytes.

`sum(result$bytes)` is the expected size of the whole object, and the
print method shows it by step.

## Details

Every block of posterior draws is a matrix of numbers with one row per
kept draw, and each number takes eight bytes. How many draws are kept is
set by `iterations`, how many columns each block has by the
specification: the number of coefficients, whether they vary over time,
the error distribution, variable selection and the rank of a VEC model.
**Time-varying parameters and stochastic volatility multiply the size of
a block by the number of periods in the sample**, which is how a model
that is small as a constant coefficient model can need gigabytes.

A few blocks depend on the priors: the inclusion indicators of the
covariance block, the draws of a non-centred prior (`omega_v`) and a
drawn autocorrelation of the cointegration space (`rho_min`). They are
counted once
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
has added the priors. Forecast draws are counted once
[`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md)
has set a horizon.

The numbers are what the object will hold in memory. R may briefly hold
a second copy while a step modifies it, so the memory a session needs
can be up to twice as large.
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
compresses what it writes, but posterior draws compress poorly, so the
size of the file is usually a little below the size in memory.

## See also

[`expected_model_size`](https://franzmohr.github.io/bvartools/reference/expected_model_size.md)
for the warnings that are based on it, and
[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
for the blocks of the object.

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# A model with constant coefficients and one with time-varying parameters
# and stochastic volatility
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 5000, burnin = 1000)
expected_model_size(model)
#> Expected size of the model once its posterior is complete: 4.8 MB 
#> 
#>                                  size   
#>  data, priors and initial values 64.3 kB
#>  add_posterior_coefficients()    1.2 MB 
#>  add_posterior_loglik()          3.6 MB 
#> 
#>  element                      draws columns size  
#>  posterior$a$coeffs           5000  21      840 kB
#>  posterior$u_sigma_inv$coeffs 5000   9      360 kB
#>  posterior$loglik             5000  89      3.6 MB

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          tvp = TRUE, error = "sv",
                          iterations = 5000, burnin = 1000)
expected_model_size(model)
#> Expected size of the model once its posterior is complete: 122.1 MB 
#> 
#>                                  size    
#>  data, priors and initial values 64.3 kB 
#>  add_posterior_coefficients()    118.4 MB
#>  add_posterior_loglik()          3.6 MB  
#> 
#>  element                      draws columns size   
#>  posterior$a$coeffs           5000  1869    74.8 MB
#>  posterior$a$sigma            5000    21    840 kB 
#>  posterior$u_sigma_inv$coeffs 5000   801    32 MB  
#>  posterior$u_sigma_inv$sigma  5000     3    120 kB 
#>  posterior$u_omega_inv$coeffs 5000   267    10.7 MB
#>  posterior$loglik             5000    89    3.6 MB 

# The total in bytes
sum(expected_model_size(model)$bytes)
#> [1] 122064264
```
