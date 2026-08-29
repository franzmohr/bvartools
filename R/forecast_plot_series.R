#' Prepare Forecast Data for Plotting
#'
#' Combines the observed values of the endogenous variables of a model with the
#' credible bands of their forecasts.
#'
#' @param x an object of class 'bvarprd'.
#' @param n_pre number of observations that precede the forecasts. If \code{NULL}
#' (default), all available observations are used.
#' @param ci interval used to calculate the credible bands of the forecasts.
#'
#' @return A named list of time-series objects. Each object contains the observed
#' series, the lower bound of the credible band, the median forecast and the
#' upper bound of the credible band of the respective endogenous variable.
#'
#' @noRd
forecast_plot_series <- function(x, n_pre = NULL, ci = 0.95) {

  y <- x[["y"]]
  tt <- nrow(y)
  var_names <- dimnames(y)[[2]]

  ci_low <- (1 - ci) / 2
  fcst <- apply(x[["fcst"]], c(1, 2), stats::quantile, probs = c(ci_low, .5, 1 - ci_low))

  ts_temp <- stats::tsp(y)
  fcst_end <- ts_temp[2] + dim(fcst)[2] / ts_temp[3]

  result <- list()
  for (i in var_names) {
    temp <- stats::ts(t(fcst[,, i]), start = ts_temp[2] + 1 / ts_temp[3], frequency = ts_temp[3])
    temp <- cbind(y[, i], temp)
    # Connect the observed series with the forecasts
    temp[tt, 2:4] <- y[tt, i]
    if (!is.null(n_pre)) {
      if (n_pre < tt) {
        temp <- temp[-c(1:(tt - n_pre)), ]
      }
      temp <- stats::ts(temp, end = fcst_end, frequency = ts_temp[3])
    }
    result[[i]] <- temp
  }

  return(result)
}
