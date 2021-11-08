#' Forecast Error Variance Decomposition
#'  
#' A generic function used to calculate forecast error varianc decompositions.
#' 
#' @param object an object of class \code{"bvar"}.
#' @param ... arguments passed forward to method.
#' 
#' @export
fevd <- function (object, ...) {
 UseMethod("fevd")
}