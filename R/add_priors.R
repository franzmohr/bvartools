#' Add Priors to Bayesian Models
#'  
#' A generic function used to generate prior specifications for a list of models. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of class \code{"bvarmodel"} or \code{"bvecmodel"}.
#' @param ... arguments passed forward to method.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Obtain data matrices
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' @export
add_priors <- function (object, ...) {
 UseMethod("add_priors")
}