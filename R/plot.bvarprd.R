#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class \code{"bvarprd"}.
#' 
#' @param x an object of class "bvarprd", usually, a result of a call to \code{\link{predict.bvar}}.
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
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Calculate forecasts
#' pred <- predict(object, new_D = rep(1, 10))
#' 
#' # Plot forecasts
#' plot(pred)
#' 
#' @export
plot.bvarprd <- function(x, ...) {
  y <- x$y
  tt <- nrow(y)
  var_names <- dimnames(y)[[2]]
  
  graphics::par(mfcol = c(length(var_names), 1))
  for (i in var_names) {
    n_ahead <- nrow(x$fcst[[i]])
    temp <- cbind(y[, i], x$fcst[[i]])
    temp[tt, 2:4] <- y[tt, i]
    stats::plot.ts(temp, plot.type = "single", lty = c(1, 2, 1, 2), main = i, ylab = "")
  }
  graphics::par(mfcol = c(1, 1))
}