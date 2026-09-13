#' Posterior Simulation of Model Coefficients
#' 
#' Generic function used for posterior simulation of model coefficients.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_posterior_coefficients.bvarmodel}},
#' \code{\link{add_posterior_coefficients.bvecmodel}},
#' \code{\link{add_posterior_coefficients.expandingwindow}},
#' \code{\link{add_posterior_coefficients.externalforecast}},
#' \code{\link{add_posterior_coefficients.modellist}}.
#'
#' @export
add_posterior_coefficients <- function(object, ...){
  UseMethod("add_posterior_coefficients")
}
