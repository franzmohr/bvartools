#' Add Forecast Input Data
#'
#' Generates and adds data matrices for forecast simulation to the elements of
#' an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with deterministic data. If not
#' specified, the function will try to identify the deterministic terms
#' automatically. If this is not successful, an error message we be returned.
#' @param exogen a time-series object with unmodelled, non-deterministic data.
#' See 'Details'.
#' @param ... arguments passed forward to method.
#'
#' @return A list of class 'bvarmodel'.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
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
#' @method add_forecast_input bvarmodel
add_forecast_input.bvarmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL, ...){
  
  fcst_input <- prepare_forecast_input(object, n_ahead = n_ahead, deterministic = deterministic, exogen = exogen, ...)
  object[["model"]][["h"]] <- fcst_input[["h"]]
  object[["data"]][["forecast"]] <- list("z" = fcst_input[["z"]])
  
  return(object)
}