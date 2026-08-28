
#' @export
#' @method get_forecast_errors bvarmodel
get_forecast_errors.bvarmodel <- function(object, ...) {
  return(object[["posterior"]][["forecast_errors"]])
}