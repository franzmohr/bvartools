#' Forecast Error Variance Decomposition
#'  
#' A generic function used to calculate forecast error variance decompositions.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
fevd <- function (object, ...) {
 UseMethod("fevd")
}