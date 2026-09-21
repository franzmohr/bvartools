#' Plotting Draws of a Bayesian Time Series Models
#' 
#' A plot function for objects of class 'modellist'.
#' 
#' @param x an object of class 'modellist'.
#' @param ... arguments passed forward to other methods.
#'
#' @return \code{x}, invisibly. The function is called for its side effect, one
#' plot per model in the list, drawn by the plot method of that model.
#'
#' @export
plot.modellist <- function(x, ...) {

  for (i in seq_along(x)) {
    plot(x[[i]], ...)
  }

  invisible(x)
}


