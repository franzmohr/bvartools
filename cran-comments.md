## Resubmission
This is a resubmission. In this version I have:

* Addressed the valgrind issue in post_coint_kls_sur from the last release (faulty example)
* Introduced checks in post_coint_kls and post_coint_kls_sur to prevent misspecification

In the previous submission, which updates an already published package, I

* updated the vignettes.
* Added function `ssvs_prior`.
* added function `minnesota_prior`.
* enhanced the functionality and error handling of the functions irf, post_coint_kls_sur, bvar, plot.bvarfevd and predict.

## Test environments
ubuntu 14.04 (on travis-ci), R-devel, R 3.6.0, R 3.5.3, R 3.4.4, R 3.3.3
win-builder (devel)

## R CMD check results

* NOTE: checking installed package size ... installed size is  9.1Mb sub-directories of 1Mb or more: libs 8.3Mb
