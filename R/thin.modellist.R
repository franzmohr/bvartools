#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws of the elements in an object of class 'modellist'.
#' 
#' @param x an object of class 'modellist'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'modellist'.
#' 
#' @export
thin.modellist <- function(x, thin = 10, ...) {
  
  for (i in 1:length(x)) {
    x[[i]] <- thin(x[[i]], thin = thin, ...)
  }
 
  return(x)
}