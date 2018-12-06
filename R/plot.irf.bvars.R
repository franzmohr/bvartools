#' Plotting Impulse Responses of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "irf.bvars".
#' 
#' @param x an object of class "irf.bvars", usually, a result of a call to \code{\link{irf}}.
#' @param ... further graphical parameters.
#' 
#' @export
plot.irf.bvars <- function(x, ...) {
  stats::plot.ts(x, plot.type = "single", col = c("red", "blue", "red"), ...)
  graphics::abline(h = 0)
}