# Export to HDF5 File

Exports the content of an object of class 'expandingwindow' to multiple
HDF5 file.

## Usage

``` r
# S3 method for class 'expandingwindow'
write_to_hdf5(object, folder, ...)
```

## Arguments

- object:

  an object of class 'expandingwindow'.

- folder:

  path to the folder where the individual models in argument `object`
  should be saved.

- ...:

  further arguments passed to or from other methods.

## Examples

``` r

# Load data
data("e1")
train <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(data = train,
                          p = 0:2,
                          deterministic = "const",
                          iterations = 10,
                          burnin = 10)
# Number of iterations and burnin should be much higher.

# Use expanding window
model <- use_expanding_window(model, start = 1982.5)

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1),
                    sigma = list(df = 3, scale = 1))

# Add initial values
model <- add_initial_values(model)

# Save model
path_to_target_directory <- tempdir()
write_to_hdf5(model, folder = path_to_target_directory)

```
