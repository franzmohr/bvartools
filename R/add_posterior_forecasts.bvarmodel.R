#' Add Forecasts
#'
#' Calculates and adds forecasts to an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel', usually, the result of a call
#' to \code{\link{add_posterior_coefficients}} and \code{\link{add_forecast_input}}.
#' @param forecast_states character, what a model with time varying coefficients or
#' stochastic volatility does with them over the forecast horizon. \code{"simulate"}
#' carries each draw's random walks forward, one step per period, so that the forecasts
#' are draws from the posterior predictive distribution of the estimated model.
#' \code{"hold"} keeps the coefficients and volatilities at their values in the last
#' sample period, which gives forecasts conditional on no further drift and narrower
#' intervals, and is what earlier versions of the package did. If \code{NULL} (default),
#' the value in \code{object$model$forecast_states} is used, and \code{"simulate"} when
#' there is none. Models with constant coefficients and volatility are unaffected.
#' @param ... arguments passed forward to method.
#'
#' @return The object in \code{object} with \code{posterior$forecast$forecasts} added, a
#' \code{\link[coda]{mcmc}} object with one row per draw and \eqn{Kh} columns, stacked by
#' period: the \eqn{K} variables of the first forecast period, then those of the second,
#' and so on. \code{posterior$forecast} is the group everything the forecast periods
#' produce hangs below, the \code{errors} of \code{\link{add_forecast_errors}} beside
#' these. \code{\link[=predict.bvarmodel]{predict}} summarises them. A
#' \code{forecast_states} that was given is stored in \code{model$forecast_states}.
#'
#' Simulating the volatility forward needs the variance of the log-volatility
#' innovations, \code{posterior$u_sigma_inv$sigma}, which \code{\link{add_posterior_coefficients}}
#' stores. A stochastic volatility model fitted with an earlier version of the package
#' lacks it and stops with an error unless \code{forecast_states = "hold"}.
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
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family posterior simulation
#' @export
add_posterior_forecasts.bvarmodel <- function(object, forecast_states = NULL, ...){

  if (!is.null(forecast_states)) {
    object[["model"]][["forecast_states"]] <- match.arg(forecast_states, c("simulate", "hold"))
  }

  algorithm <- object[["model"]][["algorithm"]]
  
  if (algorithm %in% c("VarNormalAld", "VarTvpAld")) {
    stop("A quantile regression model does not forecast: the h step ahead quantile is not the ",
         "quantile of the iterated one step ahead quantiles, so a simulated path could not be ",
         "read as a quantile of anything.")
  }
  
  if (is.null(object[["model"]][["h"]])) {
    stop("Model specification does not contain forecast horizon 'h'. Consider using function add_forecast_input().")
  }
  
  
  # Either spelling of the out-of-sample regressors will do: `x` is the compact
  # layout this package writes, `z` the SUR one an object fitted by an earlier
  # version carries, which the C++ side compacts on the way in.
  if (is.null(object[["data"]][["forecast"]][["x"]]) & is.null(object[["data"]][["forecast"]][["z"]]) &
      !is.null(object[["data"]][["train"]][["z"]]) & !object[["model"]][["structural"]]) {
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
  object[["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(object[["posterior"]][["forecast"]][["forecasts"]],
                                                                  start = mcpar_temp[1], end = mcpar_temp[2], thin = mcpar_temp[3])
  
  return(object)
}
