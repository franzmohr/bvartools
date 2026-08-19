#' Add Forecast Errors
#'  
#' Generic function used to calculate forecast errors and add them to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_forecast_errors <- function (object, test_sample, ...) {
 UseMethod("add_forecast_errors")
}