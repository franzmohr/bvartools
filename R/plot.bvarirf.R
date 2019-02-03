#' Plotting Impulse Responses of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarirf".
#' 
#' @param x an object of class "bvarirf", usually, a result of a call to \code{\link{irf}}.
#' @param ... further graphical parameters.
#' 
#' @export
#' @rdname irf
plot.bvarirf <- function(x, ...) {
  x <- cbind(0, x)
  stats::plot.ts(x, plot.type = "single", lty = c(1, 2, 1, 2), ...)
}