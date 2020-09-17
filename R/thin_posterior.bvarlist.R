#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvarlist"}.
#' 
#' @param x an object of class \code{"bvarlist"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate multiple model matrices
#' model <- gen_var(e1, p = 1:2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Thin
#' object <- thin_posterior(object)
#' 
#' @return An object of class \code{"bvarlist"}.
#' 
#' @export
thin_posterior.bvarlist <- function(x, thin = 10) {
  
  for (i in 1:length(x)) {
    x[[i]] <- thin_posterior(x[[i]], thin = thin)
  }
 
  return(x)
}