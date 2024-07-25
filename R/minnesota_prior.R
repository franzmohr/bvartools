#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a model.
#'  
#' A generic function used to generate data matrices based on the
#' Minnesota prior. The function invokes particular methods which depend on the
#' class of the first argument.
#' 
#' @param object an object of class \code{"bvarmodel"} or \code{"bvecmodel"}.
#' @param ... arguments passed forward to method.
#' 
#' @export
minnesota_prior <- function (object, ...) {
 UseMethod("minnesota_prior")
}