#' Posterior Simulation of Model Coefficients
#' 
#' Generic function used for posterior simulation of model coefficients.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_posterior_coefficients <- function(object, ...){
  UseMethod("add_posterior_coefficients")
}
