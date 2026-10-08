
# The draws of the forecast paths, or NULL where a model has not been forecast.
#
# posterior$forecast is a list, matching the group the model file keeps them in:
# 'forecasts' beside the 'errors' of add_forecast_errors() and, in time, the
# draws of the predictive density. An object fitted before the forecasts moved
# in there carries them as a matrix at posterior$forecast itself, which is
# refused rather than read -- everything downstream indexes the group, and a
# silent NULL would read as "this model was never forecast".
.forecast_draws <- function(object) {

  forecast <- object[["posterior"]][["forecast"]]
  if (is.null(forecast)) {
    return(NULL)
  }
  if (!is.list(forecast)) {
    stop("Element posterior$forecast of this object holds a matrix of draws, which is ",
         "the layout of bvartools before the forecasts moved to ",
         "posterior$forecast$forecasts. Use 'add_posterior_forecasts' on the model ",
         "again to obtain them in the current layout.")
  }

  return(forecast[["forecasts"]])
}
