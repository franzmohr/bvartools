#' Expanding Window Estimation
#' 
#' Creates objects for expanding window posterior simulation.
#' 
#' @param object an object containing model specifications and input data.
#' Usually, the output of a call to  \code{\link{create_bvarmodel}} or
#' \code{\link{create_bvecmodel}}.
#' @param ... arguments passed forward to method.
#' 
#' @export
use_expanding_window <- function(object, ...) {
  UseMethod("use_expanding_window")
}