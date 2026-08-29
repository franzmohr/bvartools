#' Prepare Forecast Input
#' 
#' Guard method for objects of class 'bvecmodel'.
#' 
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#' 
#' @export
#' @method prepare_forecast_input bvecmodel
prepare_forecast_input.bvecmodel <- function(object, ...) {
  
  stop("'prepare_forecast_input' does not work directly on a 'bvecmodel' object.\n",
       "Use 'vec_to_var()' first and then use 'prepare_forecast_input' on the\n",
       "resulting 'bvarmodel' object.",
       call. = FALSE)
  
}
