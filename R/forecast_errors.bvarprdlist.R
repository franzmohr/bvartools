#' Forecast Errors
#'
#' Calculates forecast errors for a list of Bayesian models.
#'
#' @param object an object of class 'bvarprdlist', usually, the
#' result of a call to \code{\link{predict.posteriorlist}}.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'fcsterrorlist'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1_long <- diff(log(e1)) * 100
#' e1 <- window(e1_long, end = c(1976, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 1:2, deterministic = "const",
#'                           iterations = 10, burnin = 2)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Generate forecasts
#' bvar_pred <- predict(object, n_ahead = 4)
#' 
#' # Forecast errors
#' fe <- forecast_errors(bvar_pred, test_sample = e1_long)
#' 
#' @export
forecast_errors.bvarprdlist <- function(object, test_sample, ...){
  
  object <- lapply(object, forecast_errors, test_sample = test_sample, ...)
  
  class(object) <- append("fcsterrorlist", class(object))
  
  return(object)
}