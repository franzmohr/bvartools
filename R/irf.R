#' Impulse Response Function
#'  
#' A generic function used to calculate impulse response functions.
#' 
#' @param x an object of class \code{"bvar"}.
#' @param ... arguments passed forward to method.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model data
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Obtain IR
#' ir <- irf(object, impulse = "invest", response = "cons")
#' 
#' @export
irf <- function (x, ...) {
 UseMethod("irf")
}