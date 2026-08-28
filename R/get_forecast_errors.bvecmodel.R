
#' @export
#' @method get_forecast_errors bvecmodel
get_forecast_errors.bvecmodel <- function(object, ...) {
  return(object[["posterior"]][["forecast_errors"]])
}