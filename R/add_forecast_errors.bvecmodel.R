#' Add Forecast Errors
#'
#' Calculates the forecast errors of a VEC model against a test sample, in levels.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_forecasts}}.
#' @param test_sample a time-series object of the endogenous variables, in levels, that
#' covers the forecast periods.
#' @param ... further arguments passed to \code{\link{add_forecast_errors.bvarmodel}}.
#'
#' @details The forecasts of a VEC model are of the levels, so the errors are taken against
#' the levels of \code{test_sample}, exactly as for the VAR representation
#' \code{\link{vec_to_var}} would give.
#'
#' @return The object in \code{object} with \code{posterior$forecast_errors} added, as
#' described in \code{\link{add_forecast_errors.bvarmodel}}.
#'
#' @family model comparison
#' @export
#' @method add_forecast_errors bvecmodel
add_forecast_errors.bvecmodel <- function(object, test_sample, ...){

  if (is.null(object[["posterior"]][["forecast"]])) {
    stop("Object does not contain forecasts.")
  }

  # The errors are the VAR representation's, whose data and specification are
  # the ones the forecast was made against. Only the forecasts go across.
  level <- .vec_level_form(object)
  level[["model"]][["h"]] <- object[["model"]][["h"]]
  level[["posterior"]] <- list("forecast" = object[["posterior"]][["forecast"]])
  level <- add_forecast_errors(level, test_sample = test_sample, ...)

  object[["posterior"]][["forecast_errors"]] <- level[["posterior"]][["forecast_errors"]]

  return(object)
}
