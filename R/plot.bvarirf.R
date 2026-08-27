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
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- add_posterior_coefficients(model)
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

  # Defaults that may be overridden via '...'. The zero line is drawn solid,
  # the credible bands dashed and the median solid.
  dots <- list(...)
  args <- list(plot.type = "single", lty = c(1, 2, 1, 2),
               ylab = "Response", xlab = "Period")
  args <- args[!(names(args) %in% names(dots))]

  do.call(stats::plot.ts, c(list(x), args, dots))
}