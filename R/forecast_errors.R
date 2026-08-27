#' Forecast Errors
#'  
#' A generic function used to calculate forecast errors.
#' 
#' @param object an object of class 'bvarprd', 'bvarprdlist' or
#' 'expandwindbvarprdlist'. Usually, the result of a call to
#' \code{\link{predict.bvarmodel}}, \code{\link{predict.posteriorlist}},
#' or \code{\link{predict.expandwindposteriorlist}}.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @export
forecast_errors <- function (object, test_sample, ...) {
 UseMethod("forecast_errors")
}