#' Forecast Error Variance Decomposition
#'
#' A generic function used to calculate forecast error variance decompositions.
#' 
#' @param x an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{x},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{fevd.bvarmodel}}.
#'
#' @export
fevd <- function (x, ...) {
 UseMethod("fevd")
}
