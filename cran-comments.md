## Submission

This is the major update 0.3.0 announced. The package has been reorganised
around model objects: a model is built by `create_bvarmodel()` or
`create_bvecmodel()`, given priors and starting values, and then simulated,
forecast, scored and compared through methods on that object. `NEWS.md` lists
the user-visible differences, and the section *Moving from 0.3.0 to 1.0.0* in
it names every function that changed or went.

The functions 0.3.0 deprecated still exist and still emit their transition
message; nothing that worked in 0.3.0 stops working without saying so first.

## Test environments

* GitHub Actions
  * ubuntu-latest: R-devel, R-release, R-oldrel-1
  * macos-latest: R-release
  * windows-latest: R-release
* local: Ubuntu 24.04 in Docker, R 4.6.1, `--as-cran` (the `ubuntu-latest` job
  of the workflow above, reproduced)
* local: Windows 11, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

The `--as-cran` run reports three INFO lines rather than notes: the C++
specification (`CXX17`, which `src/Makevars` sets for the vendored BayesTS
core), GNU make as a `SystemRequirements`, and the installed size.

The `--as-cran` run above reports an installed size of 87.9Mb, of which `R` is
1.4Mb and `doc` 1.2Mb; the rest is `libs`. That is compiled C++ implementing
the posterior simulators of twenty algorithms and the two discounted
estimators, built unstripped. The figure on a machine that strips debug
symbols, as CRAN's builders do, is a fraction of it.

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

<!-- TODO before submitting: confirm the version currently on CRAN. This file
     says 0.3.0, on the strength of the tag that recorded the state that went
     there; check the CRAN page rather than the tag. -->
