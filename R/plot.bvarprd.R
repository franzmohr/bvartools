#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class 'bvarprd'.
#' 
#' @param x an object of class 'bvarprd', usually, a result of a call to \code{\link{predict.bvarmodel}}.
#' @param n_pre number of plotted observations that precede the forecasts. If \code{NULL} (default),
#' all available observations will be plotted.
#' @param ci interval used to calculate the credible bands of the forecasts.
#' @param ... further graphical parameters. Arguments \code{main}, \code{ylab},
#' \code{lty} and \code{plot.type} overwrite the defaults of the function.
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
#' model <- add_posterior_coefficients(model)
#' 
#' # Calculate forecasts
#' pred <- predict(model)
#' 
#' # Plot forecasts
#' plot(pred)
#' 
#' @export
plot.bvarprd <- function(x, n_pre = NULL, ci = 0.95, ...) {
  
  dots <- list(...)

  series <- forecast_plot_series(x, n_pre = n_pre, ci = ci)

  for (i in names(series)) {
    # Defaults that may be overridden via '...'
    args <- list(plot.type = "single", lty = c(1, 2, 3, 2), main = i, ylab = "")
    args <- args[!(names(args) %in% names(dots))]

    do.call(stats::plot.ts, c(list(series[[i]]), args, dots))
  }
}