# Run BayesTS on Stored Models

Returns a function that runs the standalone BayesTS executable on model
files or on directories of them, several at a time, leaving the results
in the files.

## Usage

``` r
bayests_files(executable = NULL, library_path = NULL)
```

## Arguments

- executable:

  the path to the BayesTS executable. If `NULL` (default), option
  `bvartools.bayests_executable` or, if that is not set, environment
  variable `BAYESTS_EXECUTABLE`.

- library_path:

  a character vector of directories the executable needs on its search
  path for its runtime libraries, put in front of `PATH` while it runs.
  A packaged BayesTS carries its libraries and does not need it; an
  executable from a build tree may.

## Value

A function of paths, a command, further arguments for the executable and
the number of processes to run at once. It returns the paths invisibly
and fails with the end of BayesTS's output if any of the runs does.

## Details

[`bayests_posterior`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md)
draws one model that is in the session: it writes the model to a scratch
file, runs BayesTS on it and reads the draws back. For models that are
already stored – a folder written with
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md),
which is how a model too large to hold is worked on – that is three
copies of the draws per model and a session holding them, for no
purpose. The function returned here runs BayesTS on the stored files
themselves and leaves what it produces where it wrote it, so the draws
never pass through R.

The returned function takes paths, which are model files or directories
of them, the name of a BayesTS command, and how many to run at once:

- `"posterior"`:

  draws the posterior, and with it the pointwise log-likelihood unless
  `--no-loglik` is passed in `args`.

- `"loglik"`:

  the log-likelihood of a model that has draws.

- `"forecasts"`:

  the forecasts of a model that has draws.

- `"check"`:

  reads and validates the models without running them, which is what to
  call before an estimation that will take hours.

A model is drawn with the seed in its file, so it does not matter
whether it is run on its own, as one of a directory, or through
[`bayests_posterior`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md):
the draws are the same.

## Running several at once

The executable works through a directory one model at a time and gains
nothing from `OMP_NUM_THREADS` on models of the size a country model of
a global VAR has, so the way to use a machine is to give it several
paths and a `jobs` above one.

They are run as processes of the operating system rather than on a
cluster of R workers, because an idle R worker is not cheap: on Windows
a fresh one commits about 2 GB before it does anything and about 4 GB
once it has loaded a package of this kind, which for six workers is more
memory than the six BayesTS processes they would be waiting for. Here
one R session starts the processes, each writing its output to a log of
its own, and waits for them.

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`chain_diagnostics()`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r

if (FALSE) { # \dontrun{
run <- bayests_files(executable = "/opt/bayests/bin/bayests")

# One directory of models, drawn in place
run("models/gvar/submodels/US")

# Every country, six processes at a time
run(list.dirs("models/gvar/submodels", recursive = FALSE), jobs = 6)

# Validate first, draw without the log-likelihood
run("models/gvar/submodels/US", command = "check")
run("models/gvar/submodels/US", args = "--no-loglik")
} # }
```
