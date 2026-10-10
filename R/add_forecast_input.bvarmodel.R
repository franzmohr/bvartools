#' Add Forecast Input Data
#'
#' Generates and adds data matrices for forecast simulation to the elements of
#' an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with deterministic data. If not
#' specified, the function will try to identify the deterministic terms
#' automatically. If this is not successful, an error message we be returned.
#' @param exogen a time-series object with the unmodelled, non-deterministic variables of the
#' model. Required if the model has such variables. It has to cover the \code{s} periods before
#' the first forecast period, which are the last \code{s} periods of the estimation sample, as
#' well as the \code{n_ahead} forecast periods, because the regressors of a forecast period
#' include the lags of these variables. See \code{\link[bayests]{prepare_forecast_input}}.
#' @param y_start an optional time-series object with the endogenous variables of the model,
#' by name and in the units of the estimation sample, that the forecast starts from: the
#' forecast periods are the \code{n_ahead} periods after its last observation, and the lags of
#' the first one are its last \code{p} observations. It may end after the estimation sample, so
#' that a forecast can start from a nowcast or from data that came in after the model was
#' estimated; \code{exogen} and \code{deterministic} then have to cover the new forecast
#' periods. If \code{NULL} (default), the forecast starts after the end of the estimation
#' sample. See 'Details' in \code{\link[bayests]{prepare_forecast_input}}.
#' @param ... arguments passed forward to method.
#'
#' @return The object in \code{object} with \code{model$h} set to \code{n_ahead} and
#' \code{data$forecast$x} added, the regressors of the forecast periods in the layout of
#' \code{data$train$x}, one row per period. If \code{y_start} is given,
#' \code{model$forecast_start} is the time of the first forecast period, which
#' \code{\link[=predict.bvarmodel]{predict}} dates the forecasts by, and \code{data$forecast$y}
#' holds the start values after the end of the sample, which it shows as the history the
#' forecasts continue.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#'
#' # A model estimated up to 1978Q4 whose forecast starts from the values
#' # of 1979Q1, which it was not estimated with
#' short <- create_bvarmodel(window(e1, end = c(1978, 4)), p = 2,
#'                           deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' short <- add_priors(short,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' short <- add_forecast_input(short, n_ahead = 4,
#'                             y_start = window(e1, end = c(1979, 1)))
#' short[["model"]][["forecast_start"]]
#'
#' @family posterior simulation
#' @export
#' @method add_forecast_input bvarmodel
add_forecast_input.bvarmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL,
                                         y_start = NULL, ...){

  fcst_input <- prepare_forecast_input(object, n_ahead = n_ahead, deterministic = deterministic,
                                       exogen = exogen, y_start = y_start, ...)
  object[["model"]][["h"]] <- fcst_input[["h"]]
  object[["data"]][["forecast"]] <- list("x" = fcst_input[["x"]])
  object <- .drop_forecast_results(object)

  # Where the forecast starts, if not right after the sample. Removed otherwise,
  # so that a model given start values once forecasts from the sample again.
  object[["model"]][["forecast_start"]] <- fcst_input[["start"]]
  if (!is.null(fcst_input[["y"]])) {
    object[["data"]][["forecast"]][["y"]] <- fcst_input[["y"]]
  }

  return(object)
}
