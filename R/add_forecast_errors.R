#' Add Forecast Errors
#'
#' Generic function used to calculate forecast errors and add them to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_forecast_errors.bvarmodel}},
#' \code{\link{add_forecast_errors.expandingwindow}},
#' \code{\link{add_forecast_errors.modellist}}.
#'
#' @export
add_forecast_errors <- function (object, test_sample, ...) {
 UseMethod("add_forecast_errors")
}
