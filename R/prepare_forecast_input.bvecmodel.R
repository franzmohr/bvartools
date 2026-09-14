#' Prepare Forecast Input
#'
#' Prepares the regressors of the forecast periods of a VEC model, in levels.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... arguments passed to \code{\link{prepare_forecast_input.bvarmodel}}, such as
#' \code{n_ahead}, \code{deterministic} and \code{exogen}.
#'
#' @details A VEC model is forecast in levels, so its forecast regressors are those of its
#' VAR representation. They are prepared from the data and specification
#' \code{\link{vec_to_var}} builds; the posterior draws are not converted.
#'
#' @return A list with the forecast horizon in element \code{h} and the regressors of the
#' forecast periods in element \code{x}, as returned by
#' \code{\link{prepare_forecast_input.bvarmodel}}.
#'
#' @export
#' @method prepare_forecast_input bvecmodel
prepare_forecast_input.bvecmodel <- function(object, ...) {

  prepare_forecast_input(.vec_level_form(object), ...)

}
