#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvar"}, \code{"bvarlist"} or \code{"bvec"}.
#' 
#' @param x an object of class \code{"bvar"}, \code{"bvarlist"} or \code{"bvec"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
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
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' object <- thin_posterior(object)
#' 
#' @export
thin_posterior <- function (x, thin) {
 UseMethod("thin_posterior")
}