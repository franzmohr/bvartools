#' Spillover Index
#'
#' A generic function used to calculate the connectedness measures of Diebold
#' and Yilmaz (2012) from the posterior draws of a model.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @export
spillover <- function(object, ...) {
  UseMethod("spillover")
}
