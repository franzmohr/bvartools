# Load needed functions from different packages for fast access

#' @importFrom Matrix chol
#' @importFrom Matrix crossprod
#' @importFrom Matrix solve
#' @importFrom Matrix t

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("bvartools", libpath)
}