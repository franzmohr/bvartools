#' Plotting Draws of a Bayesian Time Series Models
#' 
#' A plot function for objects of class 'modellist'.
#' 
#' @param x an object of class 'modellist'.
#' @param ... arguments passed forward to other methods.
#' 
#' 
#' @export
plot.modellist <- function(x, ...) {
  
  for (i in 1:length(x)) {
    plot(x[[i]], ...)
  }
  
}


