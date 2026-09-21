#' Add Forecast Input Data
#'
#' Prepares the input data of the forecast periods of a VEC model.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_coefficients}}.
#' @param n_ahead an integer of the forecast horizon.
#' @param deterministic,exogen the values of the deterministic terms and of the
#' unmodelled variables in the forecast periods, as for
#' \code{\link{add_forecast_input.bvarmodel}}.
#' @param ... further arguments passed to \code{\link{prepare_forecast_input}}.
#'
#' @details A VEC model is forecast in levels: it is the same model as its VAR
#' representation, and the forecast of the differences follows from that of the
#' levels. The regressors of the forecast periods are therefore those of the VAR
#' in levels, and are assembled from the data and specification that
#' \code{\link{vec_to_var}} builds. The posterior draws are not converted.
#'
#' @return The object in \code{object} with the forecast horizon in \code{model$h}
#' and the regressors of the forecast periods, in levels, in \code{data$forecast$x}.
#'
#' @family posterior simulation
#' @export
#' @method add_forecast_input bvecmodel
add_forecast_input.bvecmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL, ...){

  fcst_input <- prepare_forecast_input(object, n_ahead = n_ahead, deterministic = deterministic,
                                       exogen = exogen, ...)
  object[["model"]][["h"]] <- fcst_input[["h"]]
  object[["data"]][["forecast"]] <- list("x" = fcst_input[["x"]])

  return(object)
}



#' Add Forecasts
#'
#' Simulates forecasts of a VEC model, in levels, from its posterior draws.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_coefficients}} and \code{\link{add_forecast_input}}.
#' @param forecast_states character, what a VEC model with time varying coefficients or
#' stochastic volatility does with them over the forecast horizon. \code{"simulate"}
#' carries each draw's random walks forward, one step per period -- the loadings and
#' short-run coefficients, the cointegration vectors through their state equation, the
#' covariance block and the log-volatilities -- and rebuilds the VAR in levels from them at
#' every period, so that the forecasts are draws from the posterior predictive distribution
#' of the estimated model. \code{"hold"} keeps them at their values in the last sample
#' period. If \code{NULL} (default), the value in \code{object$model$forecast_states} is
#' used, and \code{"simulate"} when there is none. Models with constant coefficients and
#' volatility are unaffected.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The forecasts are those of the VAR in levels the VEC model implies. For a model
#' with constant coefficients they are the forecasts \code{\link{vec_to_var}} followed by
#' \code{\link{add_posterior_forecasts.bvarmodel}} gives. For a model with time varying
#' coefficients they are not: its VAR representation has no state equation of its own and
#' holds its coefficients at the last period, while this method steps the states of the VEC
#' model and converts them to levels anew in every forecast period.
#'
#' Simulating the cointegration vectors forward depends on the units of the data. A step
#' \eqn{\eta_t} of their state equation moves the error correction term by
#' \eqn{\eta_t^{\prime} w_t}, which for series far from zero -- log levels times 100,
#' say -- is of the order of their levels rather than of their variation: the random
#' walk intercept described in section 'Prior on the cointegration space' of
#' \code{\link{cointspace_prior}}, carried forward with nothing in the data to restrain
#' it. The forecast is then mostly that drift, and its intervals widen accordingly.
#' A warning names the series in the error correction term that are more than 50
#' times their standard deviation per period away from zero at the end of the sample;
#' \code{forecast_states = "hold"} is the alternative. Series near zero, such as
#' interest rates and inflation in percent, do not draw it.
#'
#' A time varying model estimated on series that \code{\link{scale_error_correction}}
#' centred or scaled cannot be simulated forward and stops with an error unless
#' \code{forecast_states = "hold"}: the state equation of its cointegration vectors
#' belongs to the transformed series, while \code{\link{rescale_error_correction}}
#' has put the vectors back on the scale of the data and moved the means into the
#' constant of the last sample period only.
#'
#' Simulating the volatility forward needs \code{posterior$u_sigma_inv$sigma}, the variance
#' of the log-volatility innovations, which \code{\link{add_posterior_coefficients}}
#' stores. A model with stochastic volatility estimated with an earlier version of the
#' package lacks it and stops with an error unless \code{forecast_states = "hold"}.
#'
#' @return The object in \code{object} with \code{posterior$forecast$forecasts} added, a
#' \code{\link[coda]{mcmc}} object with one row per draw and \eqn{Kh} columns of forecasts
#' of the levels, stacked by period. \code{predict} summarises them. A
#' \code{forecast_states} that was given is stored in \code{model$forecast_states}.
#'
#' @family posterior simulation
#' @export
#' @method add_posterior_forecasts bvecmodel
add_posterior_forecasts.bvecmodel <- function(object, forecast_states = NULL, ...){

  if (!is.null(forecast_states)) {
    object[["model"]][["forecast_states"]] <- match.arg(forecast_states, c("simulate", "hold"))
  }

  if (is.null(object[["posterior"]])) {
    stop("Argument 'object' does not contain posterior draws. Use add_posterior_coefficients first.")
  }
  if (is.null(object[["model"]][["h"]])) {
    stop("Model specification does not contain forecast horizon 'h'. Consider using function add_forecast_input().")
  }

  algorithm <- object[["model"]][["algorithm"]]
  if (is.null(algorithm)) {
    stop("Element 'model$algorithm' is missing. Was the object produced by create_bvecmodel?")
  }

  if (.is_discount(object)) {
    return(.discount_forecasts(object))
  }

  .check_simulated_coint_states(object)

  class_of_object <- class(object)

  object <- switch(algorithm,
                   VecKlgs2010 = .VecKlgs2010Forecasts(object),
                   VecNormalGamma = .VecNormalGammaForecasts(object),
                   VecNormalStochvol = .VecNormalStochvolForecasts(object),
                   VecNormalWishart = .VecNormalWishartForecasts(object),
                   VecTvpGamma = .VecTvpGammaForecasts(object),
                   VecTvpStochvol = .VecTvpStochvolForecasts(object),
                   VecTvpWishart = .VecTvpWishartForecasts(object),
                   stop("Algorithm '", algorithm, "' not supported."))

  mcpar_temp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(object[["posterior"]][["forecast"]][["forecasts"]],
                                                                  start = mcpar_temp[1], end = mcpar_temp[2],
                                                                  thin = mcpar_temp[3])

  class(object) <- class_of_object

  return(object)
}



#' Predict Method for Objects of Class bvecmodel
#'
#' Summarises the forecasts of a VEC model, in levels.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_forecasts}}.
#' @param n_ahead number of steps ahead at which to predict. If \code{NULL} (default), every
#' period simulated by \code{\link{add_posterior_forecasts}}.
#' @param ... additional arguments passed to \code{\link{predict.bvarmodel}}.
#'
#' @details The forecasts are of the levels, and so is the history they are shown with: it
#' is recovered from the differences and the error correction term the way
#' \code{\link{vec_to_var}} recovers it.
#'
#' @return A time-series object of class \code{"bvarprd"}, as returned by
#' \code{\link{predict.bvarmodel}}.
#'
#' @family posterior simulation
#' @export
#' @method predict bvecmodel
predict.bvecmodel <- function(object, n_ahead = NULL, ...){

  if (is.null(.forecast_draws(object))) {
    stop("Missing element object$posterior$forecast$forecasts. You might want to use\n'add_forecast_input' and then 'add_posterior_forecasts'\nbefore this function.")
  }

  level <- .vec_level_form(object)
  level[["model"]][["h"]] <- object[["model"]][["h"]]
  level[["posterior"]] <- list("forecast" = object[["posterior"]][["forecast"]])

  return(stats::predict(level, n_ahead = n_ahead, ...))
}



# The VAR in levels of a VEC model's data and specification, without its draws:
# where the forecast of a VEC model takes its regressors and its history from.
# The draws stay behind -- the VEC samplers forecast from their own, and
# converting them as well would only be discarded.
.vec_level_form <- function(object) {
  object[["posterior"]] <- NULL
  vec_to_var(object)
}



# How far from zero, in multiples of their typical change per period, the series
# in the error correction term of a time varying VEC model may be before
# simulating its cointegration vectors over the forecast horizon draws a
# warning. Interest rates and inflation in percent sit near one, log levels
# times 100 in the hundreds.
.coint_level_ratio_limit <- 50

# Under forecast_states = "simulate" a time varying VEC model steps its
# cointegration vectors through their state equation in every forecast period.
# A step eta moves the error correction term by eta' w, which for series far
# from zero is of the order of their levels rather than of their variation: the
# random walk intercept of 'Prior on the cointegration space' in
# ?cointspace_prior, now carried forward with nothing in the data to restrain
# it. The forecast is then the model's, but mostly that drift. It is said once,
# as a warning, and the forecast goes ahead.
#
# A model estimated on series that scale_error_correction() centred or scaled is
# refused instead. Its state equation belongs to the transformed series, while
# rescale_error_correction() has put beta back on the scale of the data and moved
# the means into the constant of the last period only, so stepping beta from
# there simulates a different model from the one that was estimated.
.check_simulated_coint_states <- function(object) {

  if (!object[["model"]][["algorithm"]] %in% c("VecTvpWishart", "VecTvpGamma", "VecTvpStochvol")) {
    return(invisible(NULL))
  }
  r <- object[["model"]][["rank"]]
  if (is.null(r) || r == 0) {
    return(invisible(NULL))
  }
  states <- object[["model"]][["forecast_states"]]
  if (!is.null(states) && states != "simulate") {
    return(invisible(NULL))
  }

  if (isTRUE(as.logical(object[["model"]][["ect_rescaled"]]))) {
    stop("This time varying VEC model was estimated on series that ",
         "'scale_error_correction' centred or scaled, so its cointegration vectors ",
         "cannot be simulated over the forecast horizon: their state equation belongs ",
         "to the transformed series, and the constant that took up the means holds ",
         "for the last period of the sample only. Use forecast_states = \"hold\".",
         call. = FALSE)
  }

  w <- object[["data"]][["train"]][["w"]]
  n_restricted <- object[["model"]][["n_restricted"]]
  if (is.null(n_restricted)) {
    n_restricted <- 0
  }
  stochastic <- seq_len(ncol(w) - n_restricted)
  values <- matrix(as.numeric(w), nrow(w))[, stochastic, drop = FALSE]
  ratio <- apply(values, 2, function(x) {
    s <- stats::sd(diff(x))
    if (!is.finite(s) || s == 0) NA_real_ else abs(x[length(x)]) / s
  })
  far <- which(ratio > .coint_level_ratio_limit)

  if (length(far) > 0) {
    names_w <- colnames(w)[stochastic]
    if (is.null(names_w)) {
      names_w <- paste0("w", stochastic)
    }
    warning("The cointegration vectors of this time varying VEC model are simulated ",
            "over the forecast horizon, and ",
            paste0(names_w[far], " (", round(ratio[far]), ")", collapse = ", "),
            " in its error correction term ", if (length(far) == 1) "is" else "are",
            " that many times ", if (length(far) == 1) "its" else "their",
            " typical change per period away from zero. A step of the vectors moves ",
            "the term by about the level of the series, so the forecasts mostly show ",
            "that drift. forecast_states = \"hold\" keeps the states at the last ",
            "period of the sample.", call. = FALSE)
  }

  invisible(NULL)
}



#' Impulse Response Function
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since impulse responses of a VEC model are
#' obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{irf}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method irf bvecmodel
irf.bvecmodel <- function(x, ...){
  .use_vec_to_var("irf")
}



#' Forecast Error Variance Decomposition
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since variance decompositions of a VEC model
#' are obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{fevd}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method fevd bvecmodel
fevd.bvecmodel <- function(x, ...){
  .use_vec_to_var("fevd")
}



#' Spillover Index
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since spillover indices of a VEC model are
#' obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{spillover}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method spillover bvecmodel
spillover.bvecmodel <- function(object, ...){
  .use_vec_to_var("spillover")
}



# The error the guards above raise. Impulse responses, variance decompositions
# and spillovers of a VEC model go through its VAR representation, so those
# functions have no method of their own for it. Without these guards they
# stopped with R's "no applicable method", which says nothing about what to do
# instead. Forecasts are the exception: they are simulated from the VEC model's
# own draws, in levels, by the methods above.
.use_vec_to_var <- function(fun) {
  stop("'", fun, "' does not work directly on a 'bvecmodel' object.\n",
       "Use 'vec_to_var()' first and then use '", fun, "' on the\n",
       "resulting 'bvarmodel' object.",
       call. = FALSE)
}
