#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a model.
#'
#' A generic function used to generate data matrices based on the
#' Minnesota prior. The function invokes particular methods which depend on the
#' class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{minnesota_prior.bvarmodel}},
#' \code{\link{minnesota_prior.bvecmodel}},
#' \code{\link{minnesota_prior.externalforecast}},
#' \code{\link{minnesota_prior.modellist}}.
#'
#' @export
minnesota_prior <- function (object, ...) {
 UseMethod("minnesota_prior")
}
