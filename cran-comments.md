## Submission

This is 1.0.0, the major update that 0.3.0, the version on CRAN, announced.
The package has been reorganised around model objects: a model is built by
`create_bvarmodel()` or `create_bvecmodel()`, given priors and starting
values, and then simulated, forecast, scored and compared through methods on
that object. `NEWS.md` lists the user-visible differences, and the section
*Moving from 0.3.0 to 1.0.0* at its top names every function that changed or
went.

The functions 0.3.0 deprecated are removed in this version, as 0.3.0 said they
would be: for one release each of them emitted a message naming its successor
or saying it had none.

## New dependency: bayests

This version imports `bayests` (>= 0.1.0), a new package by the same
maintainer that holds what `bvartools` shares with two further packages of
Bayesian time series models: the generics of the model workflow, the model
file format, lists of models and model comparison. It was submitted to CRAN
ahead of this version and is accepted there; `bvartools` re-exports the
generics, so code written for it keeps its function names.

## Test environments

* GitHub Actions
  * ubuntu-latest: R-devel, R-release, R-oldrel-1
  * macos-latest: R-release
  * windows-latest: R-release
* local: Ubuntu 24.04 in Docker, R 4.6.1, `--as-cran` (the `ubuntu-latest` job
  of the workflow above, reproduced), with bayests 0.1.0 installed from source
* local: Windows 11, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

The `--as-cran` run reports three INFO lines rather than notes: the C++
specification (`CXX17`, which `src/Makevars` sets for the vendored BayesTS
core), GNU make as a `SystemRequirements`, and the installed size.

That size is 102.6Mb, of which `R` is 1.6Mb and `doc` 1.2Mb; the rest is
`libs`. It is compiled C++ implementing the posterior simulators of fifteen
algorithms and the two discounted estimators, built unstripped. The figure on
a machine that strips debug symbols, as CRAN's builders do, is a fraction of
it.

## Reverse dependencies

`FAVAR` is the only package on CRAN that depends on this one. It imports
`bvartools` and calls `gen_var()`, `minnesota_prior()` and `post_normal()`.

This version removes `gen_var()`, so `R CMD check` on FAVAR 0.1.3 against
`bvartools` 1.0.0 gives one note, "Missing or unexported object:
'bvartools::gen_var'", and one error: FAVAR's tests call `FAVAR::FAVAR()`, which
calls `FAVAR:::BVAR()`, which stops at `gen_var()`. `minnesota_prior()` is now a
generic with methods for the model objects this package creates, so FAVAR's call
to it would no longer dispatch either. `post_normal()` is unchanged.

<!-- TODO before submitting: state here that FAVAR's maintainer was notified,
     with the date. -->

<!-- Before submitting, in this order:
     1. bayests 0.1.0 is on CRAN; confirm on its CRAN page, then state it in
        "New dependency" above (the sentence assumes it).
     2. Remove the github::franzmohr/bayests-r entries from the three workflows
        in .github/workflows and the bayests line from the README's
        development-version install.
     3. DESCRIPTION: Version 1.0.0 and Date the day of submission, and
        CITATION.cff date-released the same day. NEWS.md, CITATION.cff and
        .zenodo.json already say 1.0.0.
     4. bgvars: lower its requirement from bvartools (>= 1.0.0.9000) to
        bvartools (>= 1.0.0) in the same step.
     5. Re-run the containerised --as-cran check (r-ci-docker) on that tree and
        update "R CMD check results" and the installed size above.
     6. Notify FAVAR's maintainer (see the TODO above). -->
