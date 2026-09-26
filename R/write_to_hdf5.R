#' Export to HDF5 File
#' 
#' A generic function used to export the content of models into HDF5 files.
#' The function invokes particular methods which depend on the class
#' of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#'
#' @details Before the method is called, the size of \code{object} is compared
#' with \code{options(bvartools.size_warning)}, and a warning is given if it is
#' larger, since that is about the disk space the files will need. See
#' \code{\link{expected_model_size}}.
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
  # Once per call from the outside, as in add_posterior_coefficients().
  if (!isTRUE(.size_check[["file"]])) {
    .check_file_size(object)
    .size_check[["file"]] <- TRUE
    on.exit(.size_check[["file"]] <- FALSE, add = TRUE)
  }
  UseMethod("write_to_hdf5")
}
