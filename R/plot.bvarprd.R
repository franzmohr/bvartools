#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class "bvarprd".
#' 
#' @param x an object of class "bvarprd", usually, a result of a call to \code{\link{predict.bvar}}.
#' @param ... further graphical parameters.
#' 
#' @export
plot.bvarprd <- function(x, ...) {
  y <- x$y
  t <- nrow(y)
  var_names <- dimnames(y)[[2]]
  
  graphics::par(mfcol = c(length(var_names), 1))
  for (i in var_names) {
    n_ahead <- nrow(x$fcst[[i]])
    temp <- cbind(y[, i], x$fcst[[i]])
    temp[t, 2:4] <- y[t, i]
    stats::plot.ts(temp, plot.type = "single", lty = c(1, 2, 1, 2), main = i, ylab = "")
  }
  graphics::par(mfcol = c(1, 1))
}