#' Plotting Forecast Error Variance Decompositions of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarfevd".
#' 
#' @param x an object of class "bvarfevd", usually, a result of a call to \code{\link{fevd}}.
#' @param ... further graphical parameters.

#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model data
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Obtain FEVD
#' vd <- fevd(object, response = "cons")
#' 
#' # Plot
#' plot(vd)
#' 
#' @export
#' @rdname fevd
plot.bvarfevd <- function(x, ...) {
  par_orig <- graphics::par("mar")
  graphics::par(mar = c(5.1, 4.1, 4.1, 7.1))
  graphics::barplot(t(x), ylab = "Percentage", xlab = "Period", names.arg = stats::time(x), ...)
  graphics::par(mar = par_orig)
  legend_names <- dimnames(x)[[2]]
  graphics::legend("left", legend = legend_names, xpd = FALSE, fill = grDevices::gray.colors(NCOL(x)), inset = 1)
}