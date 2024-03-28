#' Add Initial Values for Bayesian Models
#'  
#' A generic function used to generate initial values of posterior draws. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Obtain data matrices
#' model <- create_var_model(e1, p = 2, deterministic = 2,
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @export
add_initial_values <- function (object, ...) {
 UseMethod("add_initial_values")
}