#!/usr/bin/env Rscript

## Refresh the vendored BayesTS core.
##
##   Rscript tools/update-bayests-core.R [path/to/BayesTS]
##
## Copies the two upstream directories into the two places they live here and
## reapplies the one local patch: <armadillo> becomes "bayests/arma.h", so that
## RcppArmadillo is seen first and the samplers draw from R's RNG. Everything
## this script knows is written up in src/core/VENDORED.md; if you change one,
## change the other.

args <- commandArgs(trailingOnly = TRUE)
upstream <- if (length(args) > 0) {
  args[1]
} else {
  Sys.getenv("BAYESTS_SOURCE", "D:/workspace-cpp/BayesTS")
}

if (!dir.exists(file.path(upstream, "src", "core"))) {
  stop("no BayesTS source tree at '", upstream, "'.\n",
       "  Pass the path as an argument or set BAYESTS_SOURCE.", call. = FALSE)
}

## Upstream sources that must not be copied. Each one needs a reason, and the
## reason belongs in VENDORED.md as well. Paths are relative to the directory
## being synced, which is how sync() sees them.
skip <- c(
  ## The dynamic factor models are dfmtools', not this package's. Nothing here
  ## includes these, so they compiled into the shared object without anything
  ## able to reach them. Note this does not remove any DFM from the vendored
  ## core altogether: every Dfm*Initial, Dfm*Input and Dfm*Draws and every
  ## Input::validate() method lives in inputs.h, results.h and inputs.cpp,
  ## which every model shares and which are therefore copied whole. What goes
  ## is the samplers themselves.
  ##
  ## Every DFM upstream adds has to be listed here, and a new one that is not
  ## does not merely bloat the object: dfm_support.h is skipped, so a sampler
  ## copied without it fails to compile. Upstream has four of them.
  "dfm_normal_gamma.h",              # inst/include/bayests/
  "models/dfm_normal_gamma.cpp",     # src/core/
  "models/dfm_support.h",            # src/core/
  "dfm_normal_stochvol.h",           # inst/include/bayests/
  "models/dfm_normal_stochvol.cpp",  # src/core/
  "dfm_tvp_gamma.h",                 # inst/include/bayests/
  "models/dfm_tvp_gamma.cpp",        # src/core/
  "dfm_tvp_stochvol.h",              # inst/include/bayests/
  "models/dfm_tvp_stochvol.cpp",     # src/core/

  ## The factor augmented VAR, for the same reason and by the same mechanism.
  ## It is a factor model, so it belongs to dfmtools rather than here, and
  ## favar_support.h includes the dfm_support.h skipped just above -- copying it
  ## would not link, it would not compile. Upstream has one of them.
  "favar_normal_wishart.h",          # inst/include/bayests/
  "models/favar_normal_wishart.cpp", # src/core/
  "models/favar_support.h",          # src/core/

  ## The Kalman filter a factor model is scored by. Only the four DFMs include
  ## it, so it would sit here unreachable; unlike the samplers above it would
  ## compile, since it reaches for nothing but predictive_score.h, which is
  ## copied. It is skipped anyway rather than carried as dead weight a CRAN
  ## reviewer would find listed in inst/COPYRIGHTS with nothing including it.
  "models/factor_score.h"            # src/core/
)


## Ours, not upstream's; a refresh must leave them alone.
keep <- c("arma.h", "VENDORED.md")

patch <- function(lines) {
  lines <- lines[lines != "#define ARMA_DONT_PRINT_FAST_MATH_WARNING"]
  sub("^#include <armadillo>$", "#include \"bayests/arma.h\"", lines)
}

copy_one <- function(from, to) {
  lines <- patch(readLines(from, warn = FALSE))
  old <- if (file.exists(to)) readLines(to, warn = FALSE) else character()
  if (identical(lines, old)) return("unchanged")
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, to)
  if (length(old) == 0L) "new" else "updated"
}

sync <- function(from_dir, to_dir, pattern) {
  rel <- list.files(from_dir, pattern = pattern, recursive = TRUE)
  rel <- setdiff(rel, skip)
  status <- vapply(rel, function(f) copy_one(file.path(from_dir, f),
                                             file.path(to_dir, f)), character(1))
  for (f in rel[status != "unchanged"]) {
    cat(sprintf("  %-9s %s\n", status[[f]], file.path(to_dir, f)))
  }
  cat(sprintf("  %d file(s), %d changed\n",
              length(rel), sum(status != "unchanged")))

  ## Any source here that upstream no longer has is either ours or left over.
  ## Object files from a previous build are neither.
  here <- list.files(to_dir, pattern = pattern, recursive = TRUE)
  extra <- setdiff(here, c(rel, keep, skip))
  if (length(extra) > 0) {
    cat("  no longer upstream -- delete by hand if it is not yours:\n")
    cat(paste0("    ", file.path(to_dir, extra), collapse = "\n"), "\n")
  }
}

cat("headers  <- ", file.path(upstream, "include/bayests"), "\n", sep = "")
sync(file.path(upstream, "include", "bayests"), "inst/include/bayests", "\\.h$")

cat("sources  <- ", file.path(upstream, "src/core"), "\n", sep = "")
sync(file.path(upstream, "src", "core"), "src/core", "\\.(h|cpp)$")

cat("\nskipped on purpose: ", paste(skip, collapse = ", "), "\n", sep = "")

## The patch is the whole point, so check it took. A file that still reaches
## for <armadillo> would compile and run and quietly draw from Armadillo's RNG
## instead of R's; src/bayests_rng_guard.cpp is the second line of defence.
left <- unlist(lapply(
  c(list.files("inst/include/bayests", "\\.h$", recursive = TRUE, full.names = TRUE),
    list.files("src/core", "\\.(h|cpp)$", recursive = TRUE, full.names = TRUE)),
  function(f) if (any(grepl("^#include <armadillo>$", readLines(f, warn = FALSE)))) f))
if (length(left) > 0) {
  stop("the Armadillo patch did not take in:\n  ",
       paste(left, collapse = "\n  "),
       "\nFix the rule in this script before building.", call. = FALSE)
}

## inst/COPYRIGHTS names the vendored files one by one and claims to name all of
## them, which is what a CRAN reviewer checks it against. Nothing about copying
## a file adds it to that list, so the list drifts silently -- it had three
## files missing by the time this check was written. Compare the two here, where
## the set that just changed is still at hand.
vendored <- c(
  list.files("inst/include/bayests", "\\.h$", recursive = TRUE, full.names = TRUE),
  list.files("src/core", "\\.(h|cpp)$", recursive = TRUE, full.names = TRUE))
vendored <- sort(sub("^\\./", "", vendored))

## The list is the indented block of paths; anything else in the file is prose.
copyrights <- readLines("inst/COPYRIGHTS", warn = FALSE)
listed <- sort(trimws(grep("^\\s+(inst|src)/\\S+$", copyrights, value = TRUE)))

missing <- setdiff(vendored, listed)
stale <- setdiff(listed, vendored)
if (length(missing) > 0 || length(stale) > 0) {
  stop("inst/COPYRIGHTS does not match what is vendored.\n",
       if (length(missing)) paste0("  add:    ", paste(missing, collapse = "\n  add:    "), "\n"),
       if (length(stale)) paste0("  remove: ", paste(stale, collapse = "\n  remove: "), "\n"),
       "The file claims to list every vendored source, so it has to.",
       call. = FALSE)
}
cat("\ninst/COPYRIGHTS lists all ", length(vendored), " vendored files.\n", sep = "")
cat("\nNow rebuild. Two things a refresh can break, both of which the build\n",
    "reports rather than hides:\n",
    "  * a new core source defining a symbol this package already has in src/\n",
    "    fails to link. Delete our copy and point its callers at the core one,\n",
    "    the way kalman_durbin_koopman_2002 and stochvol_ocsn_2007 went.\n",
    "  * a new core source that does not compile is upstream's problem; add it\n",
    "    to 'skip' above with a reason until it does.\n",
    "inst/COPYRIGHTS is checked against the vendored set above, so a refresh\n",
    "that added or removed a file has already refused to finish if that list\n",
    "was not updated to match.\n",
    "src/Makevars needs no edit -- it globs core/*.cpp and core/*/*.cpp -- and\n",
    "Rcpp::compileAttributes() has nothing to do here: the core exports nothing\n",
    "to R.\n", sep = "")
