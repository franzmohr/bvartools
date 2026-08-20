#' Add Log-Likelihood
#' 
#' Generic function that calculates the posterior log-likelihoods of a model.
#'  
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_posterior_loglik <- function (object, ...) {
  UseMethod("add_posterior_loglik")
}