#' Add Forecast Errors
#'
#' Calculates forecast errors and adds them to an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#'
#' @return A list of class 'bvarmodel'.
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
#' @export
add_forecast_errors.bvarmodel <- function(object, test_sample, ...){
  
  if (is.null(object[["posterior"]][["forecast"]])) {
    stop("Object does not contain forecasts.")
  }
  k <- object[["model"]][["k"]]
  draws <- nrow(object[["posterior"]][["forecast"]])
  h <- object[["model"]][["h"]]
  
  if (k == 1) {
    tsp_test_sample <- stats::tsp(test_sample)
    test_sample <- stats::ts(as.matrix(test_sample), class = c("mts", "ts", "matrix"))
    stats::tsp(test_sample) <- tsp_test_sample
  } else {
    test_sample <- test_sample[, object[["model"]][["endogen"]]]
  }
  test_sample <- stats::na.omit(test_sample)
  
  # Determine when the forecasts start
  tsp_train <- stats::tsp(object[["data"]][["train"]][["y"]])
  forecast_starts_at <- tsp_train[2] + 1 / tsp_train[3]
  
  if (forecast_starts_at %in% stats::time(test_sample)) {
    test_sample <- stats::window(test_sample, start = forecast_starts_at)
    
    if (nrow(test_sample) < h) {
      h <- nrow(test_sample)
    }
    test_sample <- as.matrix(test_sample[1:h, ])
    
    mc_stats <- coda::mcpar(object[["posterior"]][["forecast"]])
    
    # Repeat the available test data and subtract corresponding forecasts without loop
    object[["posterior"]][["forecast_errors"]] <- coda::mcmc(t(matrix(t(test_sample), h * k, draws)) - object[["posterior"]][["forecast"]][, 1:(h * k)],
                                                             start = mc_stats[1], end = mc_stats[2], thin = mc_stats[3])
  }
  
  return(object)
}