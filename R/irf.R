#' Impulse Response Function
#'  
#' A generic function used to calculate impulse response functions.
#' 
#' @param x an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
irf <- function (x, ...) {
 UseMethod("irf")
}