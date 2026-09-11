## Purpose of this release

This is a transition release. It is the functionality of 0.2.4 with the fixes
described in NEWS.md, and it precedes a version 1.0.0 that reorganises the
package around a different set of functions.

Functions that will not exist in 1.0.0 now emit a message the first time they
are used in a session, naming their successor. The message is emitted once per
function per session, is a message rather than a warning, and can be switched
off with `options(bvartools.transition.messages = FALSE)`. The package's own
examples, tests and vignettes are unaffected by it.

No reverse dependency is affected: nothing is removed or renamed in this
release, and the registered C++ entry points in `inst/include` are unchanged, so
packages that link to bvartools continue to compile and run against it.

## Test environments

* ubuntu 24.04 (on GitHub Actions): R-devel, R-4.6.1, R-4.5.3
* macOS (on GitHub Actions): R-4.6.1
* local Windows 11, R-4.6.1
* win-builder: R-devel


## R CMD check results

0 errors | 0 warnings | 0 note
