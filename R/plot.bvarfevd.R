#' Plotting Forecast Error Variance Decompositions of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarfevd".
#' 
#' @param x an object of class "fevd.bvar", usually, a result of a call to \code{\link{irf}}.
#' @param ... further graphical parameters.
#' 
#' @export
#' @rdname fevd
plot.bvarfevd <- function(x, ...) {
  graphics::barplot(t(x), ylab = "Percentage", ...)
}