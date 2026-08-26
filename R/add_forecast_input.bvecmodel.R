#' Add Forecast Input Data
#' 
#' Guard method for objects of class 'bvecmodel'.
#' 
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @export
#' @method add_forecast_input bvecmodel
add_forecast_input.bvecmodel <- function(object, ...){
  
  stop("'add_forecast_input' does not work directly on a 'bvecmodel' object.\n",
       "Use 'vec_to_var()' first and then use 'add_forecast_input' on the\n",
       "resulting 'bvarmodel' object.",
       call. = FALSE)
}