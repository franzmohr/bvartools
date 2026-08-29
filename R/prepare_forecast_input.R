#' Prepare Forecast Input
#' 
#' Generic function that generates data matrices serving as data input for downstream
#' forecasting functions.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
prepare_forecast_input <- function (object, ...) {
  UseMethod("prepare_forecast_input")
}