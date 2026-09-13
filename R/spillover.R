#' Spillover Index
#'
#' A generic function used to calculate the connectedness measures of Diebold
#' and Yilmaz (2012) from the posterior draws of a model.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{spillover.bvarmodel}},
#' \code{\link{spillover.expandingwindow}}, \code{\link{spillover.modellist}}.
#'
#' @export
spillover <- function(object, ...) {
  UseMethod("spillover")
}
