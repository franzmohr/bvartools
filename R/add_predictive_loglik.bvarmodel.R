#' Add the Log Predictive Density of a Forecast
#'
#' Scores the forecast of an object of class 'bvarmodel' against the observations
#' its horizon realised.
#'
#' @param object an object of class 'bvarmodel', usually, the result of a call to
#' \code{\link{add_posterior_forecasts}}.
#' @param test_sample a time-series object used as test data. If \code{NULL}
#' (default), the values in \code{data$test$y} of the object are used, which is
#' what a model carries after it has been scored once and what a model read from
#' a file was written with.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The draws of the log predictive density of period \eqn{T + i} are stored in
#' \code{posterior$forecast$loglik}, one row per draw and one column per scored
#' period, beside the paths in \code{posterior$forecast$forecasts} that they
#' score. The numbers are BayesTS's own: this method fills \code{data$test$y}
#' and hands the object to the same C++ that computes the score when the
#' \code{bayests} programme is run on a model file.
#'
#' Each column conditions on the observations the periods before it realised,
#' not on the path the forecast simulated. The log of the mean of the draws of a
#' column is therefore the one step ahead predictive density given everything
#' known up to that period, and those sum over the columns to
#' \eqn{\ln p(y_{T+1}, \ldots, y_{T+h} | y_{1}, \ldots, y_{T})}, the log
#' predictive likelihood of the whole realised stretch, which
#' \code{\link{selection_criteria}} reports as \code{LPL}. Scoring against a
#' simulated history instead would give the marginal density of each period on
#' its own, which is a different quantity and one whose columns could not be
#' added up.
#'
#' With the history realised rather than simulated, the regressors of the scored
#' periods do not depend on the draw. The score is then the model's own pointwise
#' log-likelihood over those periods -- the same expression as
#' \code{\link{add_posterior_loglik}} evaluates, over a different sample -- so
#' nothing is written down twice and the two cannot drift apart. What still moves
#' with the draw is everything that follows a state equation: time varying
#' coefficients, time varying error covariances and stochastic volatilities take
#' one step of their random walk per scored period, in the order and subject to
#' the same \code{model$forecast_states} that \code{\link{add_posterior_forecasts}}
#' steps them in. Under \code{"simulate"} those steps are drawn, so the score is
#' drawn too and two runs give two answers, exactly as two runs of a forecast do.
#'
#' Fewer realised periods than the horizon is not an error. The periods that are
#' there are the ones that can be scored, and the rest of the forecast is left
#' alone. A \code{test_sample} that does not reach the forecast at all leaves the
#' object unchanged, which is what a model estimated to the end of a series looks
#' like.
#'
#' Structural models are refused: their regressors include the contemporaneous
#' observations, so the realised row is not built from its lags alone, and their
#' density carries the Jacobian of \eqn{A_0} besides. So are the asymmetric
#' Laplace algorithms of quantile estimation, which do not forecast in the first
#' place.
#'
#' @return The object in \code{object} with \code{posterior$forecast$loglik}
#' added, a \code{\link[coda]{mcmc}} object with one row per draw and one column
#' per scored period, and with the values it was scored against in
#' \code{data$test$y}.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' orig <- diff(log(e1)) * 100
#' train <- window(orig, end = c(1978, 4))
#'
#' # Create model
#' model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#'
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' # Forecast the four periods that were held back
#' model <- add_forecast_input(model, n_ahead = 4)
#' model <- add_posterior_forecasts(model)
#'
#' # Score them against what those periods realised
#' model <- add_predictive_loglik(model, test_sample = orig)
#' dim(model[["posterior"]][["forecast"]][["loglik"]])
#'
#' # The log predictive likelihood is criterion "LPL"
#' selection_criteria(model)[["LPL"]]
#'
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family model comparison
#' @export
#' @method add_predictive_loglik bvarmodel
add_predictive_loglik.bvarmodel <- function(object, test_sample = NULL, ...) {

  if (is.null(.forecast_draws(object))) {
    stop("Object does not contain forecasts. Use 'add_posterior_forecasts' first.")
  }
  .check_predictive_algorithm(object)

  realised <- if (is.null(test_sample)) {
    .realised_values(object)
  } else {
    .align_test_sample(object, test_sample)
  }
  if (is.null(realised)) {
    return(object)
  }

  # The C++ entry points return a plain list, so the class is put back by hand.
  class_of_object <- class(object)

  # What a model was scored against travels with it, so that the same model can
  # be scored again without the sample being supplied a second time, and so that
  # a model written to a file carries what /data/test/y holds.
  object[["data"]][["test"]][["y"]] <- realised

  algorithm <- object[["model"]][["algorithm"]]
  object <- switch(algorithm,
                   VarNormalGamma = .VarNormalGammaScore(object),
                   VarNormalStochvol = .VarNormalStochvolScore(object),
                   VarNormalWishart = .VarNormalWishartScore(object),
                   VarTvpGamma = .VarTvpGammaScore(object),
                   VarTvpStochvol = .VarTvpStochvolScore(object),
                   VarTvpWishart = .VarTvpWishartScore(object),
                   stop("Algorithm '", algorithm, "' cannot be scored."))

  mcpar_temp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["forecast"]][["loglik"]] <-
    coda::mcmc(object[["posterior"]][["forecast"]][["loglik"]],
               start = mcpar_temp[1], end = mcpar_temp[2], thin = mcpar_temp[3])

  class(object) <- class_of_object

  return(object)
}
