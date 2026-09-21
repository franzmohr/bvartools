#' Add the Log Predictive Density of a Forecast
#'
#' Scores the forecast of an object of class 'bvecmodel' against the levels its
#' horizon realised.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_forecasts}}.
#' @param test_sample a time-series object of the endogenous variables, in
#' levels, that covers the forecast periods. If \code{NULL} (default), the values
#' in \code{data$test$y} of the object are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The forecasts of a VEC model are of the levels, so the score is of the levels
#' too, and it is taken in the VAR representation \code{\link{vec_to_var}} would
#' give: every draw is converted to its level VAR coefficients and the realised
#' row is evaluated against the lags of the realised rows before it. That is the
#' same route \code{\link{add_forecast_errors.bvecmodel}} takes and the same one
#' the forecast itself took, so the three describe one model rather than three.
#'
#' Everything else is as \code{\link{add_predictive_loglik.bvarmodel}} describes
#' it, including what each column conditions on and how a cointegration space or
#' a volatility that moves with time is carried across the horizon. A model whose
#' error correction term was scaled or centred has to be put back with
#' \code{\link{rescale_error_correction}} first, since the score is taken against
#' the levels themselves.
#'
#' @return The object in \code{object} with \code{posterior$forecast$loglik}
#' added, a \code{\link[coda]{mcmc}} object with one row per draw and one column
#' per scored period, and with the levels it was scored against in
#' \code{data$test$y}.
#'
#' @examples
#'
#' # Load data
#' data("e6")
#' e6 <- e6 * 100
#' train <- window(e6, end = c(1997, 4))
#'
#' # Create model
#' model <- create_bvecmodel(train, p = 2, r = 1, const = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 0.0001))
#'
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' # Forecast the periods that were held back and score them
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#' model <- add_predictive_loglik(model, test_sample = e6)
#'
#' dim(model[["posterior"]][["forecast"]][["loglik"]])
#'
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family model comparison
#' @export
#' @method add_predictive_loglik bvecmodel
add_predictive_loglik.bvecmodel <- function(object, test_sample = NULL, ...) {

  if (is.null(.forecast_draws(object))) {
    stop("Object does not contain forecasts. Use 'add_posterior_forecasts' first.")
  }
  .check_predictive_algorithm(object)

  realised <- if (is.null(test_sample)) {
    .realised_values(object)
  } else {
    # The levels are the VAR representation's variables, and it is that
    # representation which knows when the forecast starts and which columns of
    # the test sample are the endogenous ones.
    level <- .vec_level_form(object)
    level[["model"]][["h"]] <- object[["model"]][["h"]]
    .align_test_sample(level, test_sample)
  }
  if (is.null(realised)) {
    return(object)
  }

  # The C++ entry points return a plain list, so the class is put back by hand.
  class_of_object <- class(object)

  object[["data"]][["test"]][["y"]] <- realised

  if (.is_discount(object)) {
    object <- .discount_score(object)
    class(object) <- class_of_object
    return(object)
  }

  .check_simulated_coint_states(object)

  algorithm <- object[["model"]][["algorithm"]]
  object <- switch(algorithm,
                   VecKlgs2010 = .VecKlgs2010Score(object),
                   VecNormalGamma = .VecNormalGammaScore(object),
                   VecNormalStochvol = .VecNormalStochvolScore(object),
                   VecNormalWishart = .VecNormalWishartScore(object),
                   VecTvpGamma = .VecTvpGammaScore(object),
                   VecTvpStochvol = .VecTvpStochvolScore(object),
                   VecTvpWishart = .VecTvpWishartScore(object),
                   stop("Algorithm '", algorithm, "' cannot be scored."))

  mcpar_temp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["forecast"]][["loglik"]] <-
    coda::mcmc(object[["posterior"]][["forecast"]][["loglik"]],
               start = mcpar_temp[1], end = mcpar_temp[2], thin = mcpar_temp[3])

  class(object) <- class_of_object

  return(object)
}
