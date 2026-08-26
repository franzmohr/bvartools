#' Add Forecast Input Data
#'
#' Generates and adds data matrices for forecast simulation to the elements of
#' an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow' containing objects that can
#' be forward to their respective `add_forecast_input_data` method.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'expandingwindow'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' orig <- diff(log(e1)) * 100
#' train <- window(orig, end = c(1982, 2))
#' 
#' 
#' # Create model
#' model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' model <- use_expanding_window(model, start = 1982.25)
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#' 
#' @export
#' @method add_forecast_input expandingwindow
add_forecast_input.expandingwindow <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- add_forecast_input(object[[i]], ...)
  }
  
  return(object)
}