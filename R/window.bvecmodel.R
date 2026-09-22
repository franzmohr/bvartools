#' Time Series Windows
#' 
#' Restricts the observations in the training sample of a model of class
#' 'bvecmodel' to the specified time window.
#' 
#' @param x an object of class 'bvecmodel'.
#' @param start the start time of the period of interest.
#' @param end the end time of the period of interest.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Posterior draws that form a path over the periods of the training
#' sample are cut to the periods that remain, so that they still refer to the
#' same observations as the data. These are the coefficients and cointegration
#' vectors of a model with time varying parameters, the draws of the error term
#' that vary by period -- as under stochastic volatility -- and the pointwise
#' log-likelihood. Draws that do not vary by period, such as those of
#' \eqn{\rho}, are left as they are.
#'
#' @return An object of class 'bvecmodel'.
#'
#' @export
#' @method window bvecmodel
window.bvecmodel <- function(x, start = NULL, end = NULL, ...) {

  k <- x[["model"]][["k"]]
  orig_time <- stats::time(x[["data"]][["train"]][["y"]])

  x[["data"]][["train"]][["y"]] <- stats::window(x[["data"]][["train"]][["y"]],
                                                      start = start, end = end, ...)
  
  if (!is.null(x[["data"]][["train"]][["w"]])) {
    x[["data"]][["train"]][["w"]] <- stats::window(x[["data"]][["train"]][["w"]],
                                                        start = start, end = end, ...) 
  }
  
  if (!is.null(x[["data"]][["train"]][["x"]])) {
    x[["data"]][["train"]][["x"]] <- stats::window(x[["data"]][["train"]][["x"]],
                                                        start = start, end = end, ...) 
  }
  
  periods <- which(orig_time %in% stats::time(x[["data"]][["train"]][["y"]]))

  if (!is.null(x[["data"]][["train"]][["z"]])) {
    pos <- rep(periods * k, each = k) - k + 1 + rep(0:(k - 1), length(periods))
    x[["data"]][["train"]][["z"]] <- x[["data"]][["train"]][["z"]][pos,]
  }

  # The helpers are the VAR method's; the cointegration vectors are a path of
  # their own, and .path_widths() knows their width from the rank.
  if (!is.null(x[["posterior"]])) {
    x[["posterior"]] <- .window_posterior(x[["posterior"]], periods, length(orig_time),
                                          .path_widths(x, k), .is_discount(x))
  }
  
  return(x)
}