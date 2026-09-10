#' Forecast Error Variance Decomposition
#'  
#' A generic function used to calculate forecast error variance decompositions.
#' 
#' @param x an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
fevd <- function (x, ...) {
 UseMethod("fevd")
}