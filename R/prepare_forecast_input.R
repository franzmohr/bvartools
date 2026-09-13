#' Prepare Forecast Input
#' 
#' Generic function that generates data matrices serving as data input for downstream
#' forecasting functions.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{prepare_forecast_input.bvarmodel}},
#' \code{\link{prepare_forecast_input.bvecmodel}}.
#'
#' @export
prepare_forecast_input <- function (object, ...) {
  UseMethod("prepare_forecast_input")
}
