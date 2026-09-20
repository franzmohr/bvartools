#' Add Forecast Errors
#'
#' Calculates forecast errors and adds them to an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param test_sample a time-series object used as test data. If \code{NULL}
#' (default), the values in \code{data$test$y} of the object are used.
#' @param ... arguments passed forward to method.
#'
#' @return The object in \code{object} with \code{posterior$forecast$errors} added, a
#' \code{\link[coda]{mcmc}} object with one row per draw and \eqn{Kh} columns in the
#' order of \code{posterior$forecast$forecasts}, which it is stored beside.
#' \code{\link{selection_criteria}} summarises them.
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
#' # Add forecast errors
#' model <- add_forecast_errors(model, test_sample = orig)
#' 
#'
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family model comparison
#' @export
add_forecast_errors.bvarmodel <- function(object, test_sample = NULL, ...){

  if (is.null(.forecast_draws(object))) {
    stop("Object does not contain forecasts.")
  }

  realised <- if (is.null(test_sample)) {
    .realised_values(object)
  } else {
    .align_test_sample(object, test_sample)
  }
  if (is.null(realised)) {
    return(object)
  }

  k <- object[["model"]][["k"]]
  h <- nrow(realised)
  draws <- nrow(.forecast_draws(object))
  mc_stats <- coda::mcpar(.forecast_draws(object))

  # What a model was scored against travels with it, so that the same model read
  # back from a file can be scored again without the sample being supplied a
  # second time.
  object[["data"]][["test"]][["y"]] <- realised

  # Repeat the available test data and subtract corresponding forecasts without loop
  object[["posterior"]][["forecast"]][["errors"]] <- coda::mcmc(t(matrix(t(realised), h * k, draws)) - .forecast_draws(object)[, 1:(h * k)],
                                                               start = mc_stats[1], end = mc_stats[2], thin = mc_stats[3])

  return(object)
}
