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
#' @examples
#'
#' data("e6")
#' model <- create_bvecmodel(e6 * 100, p = 1, r = 1, const = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' recent <- window(model, start = c(1980, 1))
#' start(recent$data$train$y)
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

  .warn_empty_dummy_variables(x)

  # The helpers are the VAR method's; the cointegration vectors are a path of
  # their own, and .path_widths() knows their width from the rank.
  if (!is.null(x[["posterior"]])) {
    x[["posterior"]] <- .window_posterior(x[["posterior"]], periods, length(orig_time),
                                          .path_widths(x, k), .is_discount(x))
  }
  x <- .window_forecast(x, periods, length(orig_time))
  
  return(x)
}