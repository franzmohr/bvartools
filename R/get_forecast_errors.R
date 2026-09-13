#' Get Forecast Errors
#'
#' A generic function used to obtain the forecast errors from a model.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @export
get_forecast_errors <- function (object, ...) {
 UseMethod("get_forecast_errors")
}
