#' Forecast Error Variance Decomposition
#'  
#' A generic function used to calculate forecast error varianc decompositions.
#' 
#' @param object an object of class \code{"bvar"}.
#' @param ... arguments passed forward to method.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate models
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Obtain FEVD
#' vd <- fevd(object, response = "cons")
#' 
#' # Plot FEVD
#' plot(vd)
#' 
#' @export
fevd <- function (object, ...) {
 UseMethod("fevd")
}