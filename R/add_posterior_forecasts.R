#' Add Forecasts
#' 
#' Generic function that calculates and adds forecasts.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_posterior_forecasts <- function (object, ...) {
  UseMethod("add_posterior_forecasts")
}