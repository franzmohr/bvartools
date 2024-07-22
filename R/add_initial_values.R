#' Add Initial Values of an MCMC Chain
#' 
#' A generic function used to generate initial values of an MCMC chain. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_initial_values <- function (object, ...) {
 UseMethod("add_initial_values")
}