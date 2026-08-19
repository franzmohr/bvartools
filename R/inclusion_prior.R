#' Prior Inclusion Probabilities
#' 
#' A generic function used to generate prior inclusion probabilities as
#' required for stochastic search variable selection (SSVS) à la George et
#' al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' The function invokes particular methods which depend on the class of the
#' first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
inclusion_prior <- function(object, ...) {
  UseMethod("inclusion_prior")
}