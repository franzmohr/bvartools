#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws of the elements in an object of class 'expandingwindow'.
#' 
#' @param x an object of class 'expandingwindow'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'expandingwindow'.
#' 
#' @export
thin.expandingwindow <- function(x, thin = 10, ...) {
  
  for (i in 1:length(x)) {
    x[[i]] <- thin(x[[i]], thin = thin, ...)
  }
 
  return(x)
}