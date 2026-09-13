#' Add Forecast Input Data
#' 
#' Generic function used to generate and add data matrices for forecast simulation
#' to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_forecast_input.bvarmodel}},
#' \code{\link{add_forecast_input.bvecmodel}},
#' \code{\link{add_forecast_input.expandingwindow}},
#' \code{\link{add_forecast_input.externalforecast}},
#' \code{\link{add_forecast_input.modellist}}.
#'
#' @export
add_forecast_input <- function (object, ...) {
  UseMethod("add_forecast_input")
}
