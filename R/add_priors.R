#' Add Priors to Bayesian Models
#'  
#' A generic function used to generate prior vectors and matrices. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_priors <- function (object, ...) {
 UseMethod("add_priors")
}