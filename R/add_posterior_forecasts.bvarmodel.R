#' Add Forecasts
#'
#' Calculates and adds forecasts to an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel', usually, the result of a call
#' to \code{\link{add_posterior_coefficients}} and \code{\link{add_forecast_input}}.
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
#'
#' @export
add_posterior_forecasts.bvarmodel <- function(object, ...){
  
  algorithm <- object[["model"]][["algorithm"]]
  
  if (is.null(object[["model"]][["h"]])) {
    stop("Model specification does not contain forecast horizon 'h'. Consider using function add_forecast_input().")
  }
  
  
  if (is.null(object[["data"]][["forecast"]][["z"]]) & !is.null(object[["data"]][["train"]][["z"]]) & !object[["model"]][["structural"]]) {
    stop("Model specification does not contain input data. Consider using function add_forecast_input().")
  }
  
  if (algorithm %in% c("VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                       "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")) {
    object <- switch(algorithm,
                     VarNormalGamma = .VarNormalGammaForecasts(object),
                     VarNormalStochvol = .VarNormalStochvolForecasts(object),
                     VarNormalWishart = .VarNormalWishartForecasts(object),
                     VarTvpGamma = .VarTvpGammaForecasts(object),
                     VarTvpStochvol = .VarTvpStochvolForecasts(object),
                     VarTvpWishart = .VarTvpWishartForecasts(object))
  } else {
    stop("Algorithm not implemented yet.")
  }
  
  mcpar_temp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["forecast"]] <- coda::mcmc(object[["posterior"]][["forecast"]], start = mcpar_temp[1], end = mcpar_temp[2], thin = mcpar_temp[3])
  
  return(object)
}