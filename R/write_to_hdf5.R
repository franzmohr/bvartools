#' Export to HDF5 File
#' 
#' A generic function used to export the content of models into HDF5 files.
#' The function invokes particular methods which depend on the class
#' of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{write_to_hdf5.bvarmodel}},
#' \code{\link{write_to_hdf5.bvecmodel}},
#' \code{\link{write_to_hdf5.expandingwindow}},
#' \code{\link{write_to_hdf5.modellist}}.
#'
#' @export
write_to_hdf5 <- function (object, ...) {
  UseMethod("write_to_hdf5")
}
