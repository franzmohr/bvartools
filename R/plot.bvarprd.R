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
  t <- ncol(y)
  var_names <- dimnames(y)[[1]]
  
  graphics::par(mfcol = c(length(var_names), 1))
  for (i in var_names) {
    n_ahead <- nrow(x$fcst[[i]])
    temp <- matrix(NA, t + n_ahead, 4)
    temp[1:t, 1] <- y[i,]
    temp[t, 2:4] <- y[i, t]
    temp[(t + 1):(t + n_ahead), 2:4] <- x$fcst[[i]]
    stats::plot.ts(temp, plot.type = "single", lty = c(1, 2, 1, 2), main = i, ylab = "", ...)
  }
  graphics::par(mfcol = c(1, 1))
}