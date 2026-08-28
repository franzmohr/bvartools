#' Predict Method for Objects of Class expandingwindow
#' 
#' Forecasting Bayesian VAR objects in a list of class 'expandingwindow'.
#' 
#' @param object an object of class 'expandingwindow'.
#' @param ... additional arguments.
#' 
#' @return A time-series object of class 'expandwindbvarprdlist'.
#' 
#' @examples
#' 
#' data("us_macrodata")
#' 
#' model <- create_bvarmodel(data = us_macrodata,
#'                           p = 1,
#'                           deterministic = "none",
#'                           error = "gamma",
#'                           iterations = 10,
#'                           burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Obtain objects for expanding window estimation
#' model <- use_expanding_window(model, start = 2007)
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(shape = 3, rate = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#' 
#' model <- add_forecast_input(model, n_ahead = 10)
#' model <- add_posterior_forecasts(model)
#' 
#' # Predict
#' prd <- predict(model, n_ahead = 4)
#' 
#' @export
#' @method predict expandingwindow
predict.expandingwindow <- function(object, ...) {
  
  result <- lapply(object, predict, ...)
  
  # If output contains list elements, which are NULL, those are dropped
  pos_null <- which(unlist(lapply(result, function(x) {is.null(x)})))
  check_null <- length(pos_null) > 0
  warn <- check_null
  warn_pos <- pos_null
  while (check_null) {
    result[[pos_null[1]]] <- NULL
    pos_null <- which(unlist(lapply(result, function(x) {is.null(x)})))
    check_null <- length(pos_null) > 0
  }
  
  class(result) <- c("expandwindbvarprdlist", class(result))
  return(result)
}
