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
#' The forecast of a VEC model starts right after its estimation sample, and argument
#' \code{y_start} of the 'bvarmodel' method is refused: the VEC samplers forecast from their
#' own draws and the end of the sample, and a time varying model steps its states from
#' there. To start a forecast from values after the sample, convert the model with
#' \code{\link{vec_to_var}} and give \code{y_start} to the resulting 'bvarmodel'.
#'
#' @return A list with the forecast horizon in element \code{h} and the regressors of the
#' forecast periods in element \code{x}, as returned by
#' \code{\link{prepare_forecast_input.bvarmodel}}.
#'
#' @export
#' @method prepare_forecast_input bvecmodel
prepare_forecast_input.bvecmodel <- function(object, ...) {

  if (!is.null(list(...)[["y_start"]])) {
    stop("Argument 'y_start' is not available for a 'bvecmodel': its forecast starts right ",
         "after the estimation sample. Convert the model with vec_to_var() and give the start ",
         "values, in levels, to add_forecast_input() of the resulting 'bvarmodel'.",
         call. = FALSE)
  }

  prepare_forecast_input(.vec_level_form(object), ...)

}
