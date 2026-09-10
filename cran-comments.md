## Submission

This is a major update. The version currently on CRAN is 0.2.4. The package has
been reorganised around model objects and some results have changed; `NEWS.md`
lists the user-visible differences.

## Test environments

* GitHub Actions
  * ubuntu-latest: R-devel, R-release, R-oldrel-1
  * macos-latest: R-release
  * windows-latest: R-release
* local: Windows 11, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 1 note

The note is "Skipping checking math rendering: package 'V8' unavailable". It
describes the local machine, where V8 is not installed, rather than the package.

`R CMD check --as-cran` reports an installed size of 5.8Mb, 3.7Mb of it in
`libs`. That is compiled C++ implementing the posterior simulators.

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
