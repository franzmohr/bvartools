#' Add Forecast Input Data
#' 
#' Generic function used to generate and add data matrices for forecast simulation
#' to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_forecast_input_data <- function (object, ...) {
  UseMethod("add_forecast_input_data")
}