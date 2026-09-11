# Pre-compile the expensive vignettes.
#
# The .Rmd files in this folder are *generated*: they contain the results of the
# code rather than the code itself, so that `R CMD build` has nothing left to
# evaluate. The real sources are the matching .Rmd.orig files, and this script
# turns the latter into the former. Both the generated .Rmd files and the
# figures under vignettes/figures/ are committed.
#
# Run this after changing a .orig source, or after a change to the package that
# alters vignette output:
#
#   Rscript vignettes/precompile.R                 # all vignettes
#   Rscript vignettes/precompile.R quantile-var    # just one
#
# The chunks call `library(bvartools)`, so what gets baked into the output is
# the *installed* package, not the sources in this working tree. Install first:
#
#   R CMD INSTALL .
#
# Posterior simulation takes a while, which is the whole reason this script
# exists. Neither this file nor the .orig sources are shipped in the built
# package -- see .Rbuildignore.

precompile <- function(name) {
  wd <- setwd("vignettes")
  on.exit(setwd(wd), add = TRUE)

  input <- paste0(name, ".Rmd.orig")
  if (!file.exists(input)) {
    stop("no such vignette source: vignettes/", input, call. = FALSE)
  }

  # A chunk that fails must fail the build. knitr's default is to record the
  # error in the output and carry on, which quietly ships a vignette whose
  # every step after the first mistake is an error message. Chunks that mean
  # to show an error still set error = TRUE themselves.
  knitr::opts_chunk$set(error = FALSE)

  # Every vignette writes to figures/ under its own prefix, so dropping that
  # prefix leaves no stale plots behind when chunks are renamed or removed.
  dir.create("figures", showWarnings = FALSE)
  unlink(list.files("figures", pattern = paste0("^", name, "-"), full.names = TRUE))

  message("Knitting ", input, " ...")
  started <- Sys.time()
  knitr::knit(input, paste0(name, ".Rmd"))
  message(
    "Done in ",
    format(round(difftime(Sys.time(), started), 1)),
    "\n"
  )
}

vignettes <- commandArgs(trailingOnly = TRUE)
if (length(vignettes) == 0) {
  vignettes <- c(
    "bvartools",
    "bvec",
    "minnesota-prior",
    "ssvs",
    "quantile-var",
    "tvp-sv-var",
    "tvp-sv-vec",
    "model-comparison",
    "horse-races"
  )
}

for (vignette in vignettes) {
  precompile(vignette)
}
