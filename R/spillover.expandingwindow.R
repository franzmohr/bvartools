#' Rolling Spillover Index
#'
#' Produces the total spillover index of Diebold and Yilmaz (2012) over the
#' windows of an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow' whose models carry posterior
#' draws.
#' @param ... arguments passed forward to \code{\link{spillover.bvarmodel}}.
#'
#' @details The index over a moving sample is the chart this literature is built
#' around, and \code{\link{use_expanding_window}} already produces the models it
#' needs. Each window contributes one value, dated at its own last observation,
#' so the result reads as the connectedness of the system as it was known up to
#' that date.
#'
#' Note that the windows expand rather than roll: every one of them starts at the
#' beginning of the sample, so later values are estimated on more data. Diebold
#' and Yilmaz use a fixed width instead, which trades that away for a constant
#' amount of information per point.
#'
#' @return A time-series object of class 'bvarspilloverts' with the median of the
#' total index and the bounds of its credible interval, one row per window.
#' Windows whose estimation failed are \code{NA}.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#'
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 1, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # Estimate over expanding windows
#' windows <- use_expanding_window(model, start = c(1976, 1))
#' windows <- lapply(windows, function(w) {
#'   add_posterior_coefficients(add_initial_values(w))
#' })
#' class(windows) <- append("expandingwindow", class(windows))
#'
#' # The index as the sample grows
#' si <- spillover(windows, n_ahead = 5)
#' plot(si)
#'
#' @references
#'
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive
#' directional measurement of volatility spillovers. \emph{International Journal of
#' Forecasting, 28}(1), 57--66.
#'
#' @export
spillover.expandingwindow <- function(object, ...) {

  n_window <- length(object)
  if (n_window == 0) {
    stop("Argument 'object' contains no models.")
  }

  result <- matrix(NA_real_, n_window, 3)
  dates <- rep(NA_real_, n_window)

  for (i in seq_len(n_window)) {
    window_i <- object[[i]]

    # Each window is dated at its own last observation: that is the information
    # set the value belongs to.
    y <- window_i[["data"]][["train"]][["y"]]
    if (stats::is.ts(y)) {
      dates[i] <- max(stats::time(y))
    } else {
      dates[i] <- i
    }

    if (isTRUE(window_i[["error"]])) {
      next
    }
    out <- try(spillover(window_i, ...), silent = TRUE)
    if (inherits(out, "try-error")) {
      next
    }

    total <- out[["total"]]
    if (is.matrix(total) && !is.null(rownames(total))) {
      result[i, ] <- c(total[1, 1],
                       total[rownames(total) == "50%", 1],
                       total[nrow(total), 1])
    } else {
      ci <- out[["specification"]][["ci"]]
      result[i, ] <- stats::quantile(total, probs = c((1 - ci) / 2, .5,
                                                      1 - (1 - ci) / 2))
    }
  }

  colnames(result) <- c("lower", "median", "upper")

  frequency <- 1
  y <- object[[1]][["data"]][["train"]][["y"]]
  if (stats::is.ts(y)) {
    frequency <- stats::frequency(y)
  }
  result <- stats::ts(result, start = dates[1], frequency = frequency)

  class(result) <- append("bvarspilloverts", class(result))

  return(result)
}

#' Plot a Rolling Spillover Index
#'
#' @param x an object of class 'bvarspilloverts'.
#' @param ... further arguments passed to \code{\link[stats]{plot.ts}}.
#'
#' @return \code{NULL}, invisibly. Called for its side effect.
#'
#' @export
plot.bvarspilloverts <- function(x, ...) {

  values <- unclass(x)
  attr(values, "tsp") <- NULL

  dots <- list(...)
  args <- list(ylab = "Percentage", xlab = "",
               main = "Total spillover index",
               plot.type = "single", lty = c(2, 1, 2),
               col = c("grey50", "black", "grey50"))
  args <- args[!(names(args) %in% names(dots))]

  do.call(stats::plot.ts, c(list(stats::ts(values, start = stats::start(x),
                                           frequency = stats::frequency(x))),
                            args, dots))

  invisible(NULL)
}
