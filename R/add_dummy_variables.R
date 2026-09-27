#' Dummy Variables
#'
#' A generic function used to add dummy variables to the deterministic terms of
#' a model.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_dummy_variables.bvarmodel}},
#' \code{\link{add_dummy_variables.bvecmodel}},
#' \code{\link{add_dummy_variables.expandingwindow}},
#' \code{\link{add_dummy_variables.modellist}}.
#'
#' @export
add_dummy_variables <- function(object, ...) {
  UseMethod("add_dummy_variables")
}
