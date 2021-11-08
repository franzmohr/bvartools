#' Impulse Response Function
#'  
#' A generic function used to calculate impulse response functions.
#' 
#' @param x an object of class \code{"bvar"}.
#' @param ... arguments passed forward to method.
#' 
#' @export
irf <- function (x, ...) {
 UseMethod("irf")
}