#' Plotting Impulse Responses of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarirf".
#' 
#' @param x an object of class "bvarirf", usually, a result of a call to \code{\link{irf}}.
#' @param ... further graphical parameters.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model data
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Add posterior specifications
#' object <- draw_posterior(model)
#' 
#' # Calculate IR
#' ir <- irf(object, impulse = "invest", response = "cons")
#' 
#' # Plot IR
#' plot(ir)
#' 
#' @export
#' @rdname irf
plot.bvarirf <- function(x, ...) {
  if (ncol(x) != 3) {
    stop("Cannot handle output of function 'irf' when keep_draws = TRUE.")
  }
  x <- cbind(0, x)
  stats::plot.ts(x, plot.type = "single", lty = c(1, 2, 1, 2), ...)
}