#' Add Forecasts
#'
#' Calculates and adds forecasts to an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel', usually, the result of a call
#' to \code{\link[bayests]{add_posterior_coefficients}} and \code{\link[bayests]{add_forecast_input}}.
#' @param forecast_states character, what a model with time varying coefficients or
#' stochastic volatility does with them over the forecast horizon. \code{"simulate"}
#' carries each draw's random walks forward, one step per period, so that the forecasts
#' are draws from the posterior predictive distribution of the estimated model.
#' \code{"hold"} keeps the coefficients and volatilities at their values in the last
#' sample period, which gives forecasts conditional on no further drift and narrower
#' intervals, and is what earlier versions of the package did. If \code{NULL} (default),
#' the value in \code{object$model$forecast_states} is used, and \code{"simulate"} when
#' there is none. Models with constant coefficients and volatility are unaffected.
#' @param scenario an optional data frame with columns \code{period}, \code{variable} and
#' \code{value}, one row per value a variable is held at in a forecast period: a
#' conditional forecast. \code{period} counts forecast periods from one and \code{variable}
#' is a name or a position among the endogenous variables. It is stored in
#' \code{data$forecast$constraints}, which is where the forecast reads it from; set that
#' element to \code{NULL} to forecast without it again. Available for models with constant
#' coefficients and \code{error = "wishart"} or \code{"gamma"}, and for a grid of quantiles.
#' @param forecast_quantile for a grid of quantiles (see \code{quantile_grid} in
#' \code{\link{create_bvarmodel}}), a level in \eqn{(0, 1)} every draw of the forecast is taken
#' at, which gives quantile paths rather than draws from the predictive distribution: the
#' paths of Chavleishvili and Manganelli (2019), whose difference with and without a
#' \code{scenario} is an impulse response at that quantile. If \code{NULL} (default), the
#' value in \code{object$model$forecast_quantile} is used, and the levels are drawn at random
#' when there is none. Zero removes a stored one.
#' @param ... arguments passed forward to method.
#'
#' @return The object in \code{object} with \code{posterior$forecast$forecasts} added, a
#' \code{\link[coda]{mcmc}} object with one row per draw and \eqn{Kh} columns, stacked by
#' period: the \eqn{K} variables of the first forecast period, then those of the second,
#' and so on. \code{posterior$forecast} is the group everything the forecast periods
#' produce hangs below, the \code{errors} of \code{\link[bayests]{add_forecast_errors}} beside
#' these. \code{\link[=predict.bvarmodel]{predict}} summarises them. A
#' \code{forecast_states} that was given is stored in \code{model$forecast_states}.
#'
#' Simulating the volatility forward needs the variance of the log-volatility
#' innovations, \code{posterior$u_sigma_inv$sigma}, which \code{\link[bayests]{add_posterior_coefficients}}
#' stores. A stochastic volatility model fitted with an earlier version of the package
#' lacks it and stops with an error unless \code{forecast_states = "hold"}.
#'
#' @details A forecast that starts after the end of the estimation sample, from argument
#' \code{y_start} of \code{\link[bayests]{add_forecast_input}}, takes its lags from the start values
#' and its coefficients from the sample. \strong{Three things that start from the last period
#' of the sample are refused for it}: a \code{scenario}, which is drawn from the end of the
#' sample; \code{forecast_states = "simulate"} for a model with time varying coefficients or
#' stochastic volatility, which would not step the states through the periods up to the start
#' values (\code{"hold"} is available); and the discounted model, for the same reason.
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
#' @references
#'
#' Chavleishvili, S., & Manganelli, S. (2019). Forecasting and stress testing with quantile
#' vector autoregression. \emph{ECB Working Paper}, 2330.
#'
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family posterior simulation
#' @export
add_posterior_forecasts.bvarmodel <- function(object, forecast_states = NULL, scenario = NULL,
                                              forecast_quantile = NULL, ...){

  if (!is.null(forecast_states)) {
    object[["model"]][["forecast_states"]] <- match.arg(forecast_states, c("simulate", "hold"))
  }

  algorithm <- object[["model"]][["algorithm"]]
  grid <- !is.null(object[["model"]][["quantiles"]])
  
  if (algorithm %in% c("VarNormalAld", "VarTvpAld") && !grid) {
    stop("A quantile regression model does not forecast: the h step ahead quantile is not the ",
         "quantile of the iterated one step ahead quantiles, so a simulated path could not be ",
         "read as a quantile of anything. A grid of quantiles does; see 'quantile_grid' in ",
         "create_bvarmodel().")
  }

  if (!is.null(forecast_quantile)) {
    if (!grid) {
      stop("Argument 'forecast_quantile' needs a grid of quantiles; see 'quantile_grid' in ",
           "create_bvarmodel().")
    }
    if (!is.numeric(forecast_quantile) || length(forecast_quantile) != 1 ||
        is.na(forecast_quantile) || forecast_quantile < 0 || forecast_quantile >= 1) {
      stop("Argument 'forecast_quantile' must be a single value in (0, 1), or zero.")
    }
    object[["model"]][["forecast_quantile"]] <- if (forecast_quantile == 0) NULL else forecast_quantile
  }

  if (!is.null(scenario)) {
    object[["data"]][["forecast"]][["constraints"]] <-
      .scenario_constraints(scenario, object[["model"]][["endogen"]], object[["model"]][["h"]])
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

  .check_shifted_forecast(object)

  # Drawn i.i.d. from the closed form rather than carried along a chain, so the
  # paths come back unlabelled and return here.
  if (.is_discount(object)) {
    return(.discount_forecasts(object))
  }

  if (algorithm %in% c("VarNormalAld", "VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                       "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")) {
    object <- switch(algorithm,
                     VarNormalAld = .VarNormalAldForecasts(object),
                     VarNormalGamma = .VarNormalGammaForecasts(object),
                     VarNormalStochvol = .VarNormalStochvolForecasts(object),
                     VarNormalWishart = .VarNormalWishartForecasts(object),
                     VarTvpGamma = .VarTvpGammaForecasts(object),
                     VarTvpStochvol = .VarTvpStochvolForecasts(object),
                     VarTvpWishart = .VarTvpWishartForecasts(object))
  } else {
    stop("Algorithm not implemented yet.")
  }
  
  mcpar_temp <- .draws_mcpar(object)
  object[["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(object[["posterior"]][["forecast"]][["forecasts"]],
                                                                  start = mcpar_temp[1], end = mcpar_temp[2], thin = mcpar_temp[3])
  
  return(object)
}

# A forecast from start values after the end of the sample (argument 'y_start'
# of add_forecast_input()) takes its lags from them and everything else from the
# sample. Where the samplers carry something else on from the last period of the
# sample, the forecast would silently start there instead, so it is refused.
.check_shifted_forecast <- function(object) {

  if (.forecast_gap(object) <= 0) {
    return(invisible(NULL))
  }

  if (!is.null(object[["data"]][["forecast"]][["constraints"]])) {
    .refuse_shifted_forecast(object, "a scenario",
                             paste0("a conditional forecast is drawn from the end of the estimation ",
                                    "sample, not from the start values. Add the periods of the start ",
                                    "values to the data the model is estimated with, or condition on ",
                                    "unmodelled variables through argument 'exogen'."))
  }

  if (.is_discount(object)) {
    .refuse_shifted_forecast(object, "a forecast of the discounted model",
                             paste0("its coefficients drift from the last period of the sample, one ",
                                    "step per period, and the periods up to the start values are ",
                                    "not counted. Estimate the model on a sample that includes them."))
  }

  algorithm <- object[["model"]][["algorithm"]]
  states <- object[["model"]][["forecast_states"]]
  if (is.null(states)) {
    states <- "simulate"
  }
  if (grepl("Tvp|Stochvol", algorithm) && states == "simulate") {
    .refuse_shifted_forecast(object, "forecast_states = \"simulate\"",
                             paste0("the time varying coefficients or volatilities would be carried ",
                                    "forward from the last period of the sample as if the forecast ",
                                    "started there, without the steps of the periods up to the start ",
                                    "values. Use forecast_states = \"hold\", which keeps them at the ",
                                    "last period of the sample, or estimate the model on a sample ",
                                    "that includes the start values."))
  }

  invisible(NULL)
}

# The pins of a conditional forecast, as the constraint set the core reads from
# data$forecast$constraints: one row, one entry and weight one per pinned value,
# all in the hard group 0.
.scenario_constraints <- function(scenario, endogen, h) {
  if (!is.data.frame(scenario) || !all(c("period", "variable", "value") %in% names(scenario))) {
    stop("Argument 'scenario' must be a data frame with columns 'period', 'variable' and 'value'.")
  }
  n <- nrow(scenario)
  if (n == 0) {
    stop("Argument 'scenario' has no rows.")
  }
  variable <- scenario[["variable"]]
  if (is.factor(variable)) {
    variable <- as.character(variable)
  }
  if (is.character(variable)) {
    pos <- match(variable, endogen)
    if (anyNA(pos)) {
      stop("Argument 'scenario' names variables the model does not have: ",
           paste(unique(variable[is.na(pos)]), collapse = ", "), ".")
    }
    variable <- pos
  }
  period <- scenario[["period"]]
  value <- as.numeric(scenario[["value"]])
  if (!is.numeric(variable) || anyNA(variable) || any(variable != round(variable)) ||
      any(variable < 1 | variable > length(endogen))) {
    stop("Column 'variable' of argument 'scenario' must name or number endogenous variables.")
  }
  if (is.null(h)) {
    stop("Model specification does not contain forecast horizon 'h'. Consider using function add_forecast_input().")
  }
  if (!is.numeric(period) || anyNA(period) || any(period != round(period)) ||
      any(period < 1 | period > h)) {
    stop("Column 'period' of argument 'scenario' must hold forecast periods from 1 to ", h, ".")
  }
  if (anyNA(value) || any(!is.finite(value))) {
    stop("Column 'value' of argument 'scenario' must be finite.")
  }
  if (anyDuplicated(data.frame(period, variable)) > 0) {
    stop("Argument 'scenario' pins a variable twice in the same period.")
  }
  list(value = value,
       group = rep(0, n),
       row = seq_len(n),
       period = as.numeric(period),
       variable = as.numeric(variable),
       weight = rep(1, n))
}
