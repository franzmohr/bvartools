#' Add Forecasts
#'
#' Calculates and adds forecasts to the elements of an object of class 'modellist'.
#'
#' @param object an object of class 'modellist', usually, the result of a call
#' to \code{\link{add_posterior_coefficients}} and \code{\link{add_forecast_input}}.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'modellist'.
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
#' model <- create_bvarmodel(data = train, p = 0:2, deterministic = "const",
#'                           iterations = 10, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#' 
#' # Add forecasts
#' model <- add_posterior_forecasts(model)
#' 
#' @export
add_posterior_forecasts.modellist <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- add_posterior_forecasts(object[[i]], ...)
  }
  
  return(object)
}