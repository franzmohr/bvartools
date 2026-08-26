#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class 'bvarprd'.
#' 
#' @param x an object of class 'bvarprd', usually, a result of a call to \code{\link{predict.bvarmodel}}.
#' @param n_pre number of plotted observations that precede the forecasts. If \code{NULL} (default),
#' all available observations will be plotted.
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

  y <- x[["y"]]
  k <- ncol(y)
  tt <- nrow(y)
  var_names <- dimnames(y)[[2]]
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(x[["fcst"]], c(1, 2), stats::quantile, probs = c(ci_low, .5, ci_high))
  result <- c()
  for (i in 1:k) {
    result <- c(result, list(stats::ts(t(temp[,, i]))))
  }
  names(result) <- dimnames(y)[[2]]
  
  ts_temp <- stats::tsp(y)
  for (i in 1:k) {
    result[[i]] <- stats::ts(result[[i]], start = ts_temp[2] + 1 / ts_temp[3], frequency = ts_temp[3])
  }
  
  for (i in var_names) {
    n_ahead <- nrow(result[[i]])
    temp <- cbind(y[, i], result[[i]])
    temp[tt, 2:4] <- y[tt, i]
    if (!is.null(n_pre)) {
      if (n_pre < tt) {
        temp <- temp[-c(1:(tt - n_pre)), ]
      }
      temp <- stats::ts(temp, end = stats::tsp(result[[i]])[2], frequency = stats::tsp(result[[i]])[3])
    }
    # Defaults that may be overridden via '...'
    args <- list(plot.type = "single", lty = c(1, 2, 3, 2), main = i, ylab = "")
    args <- args[!(names(args) %in% names(dots))]

    do.call(stats::plot.ts, c(list(temp), args, dots))
  }
}