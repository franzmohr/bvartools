#' Predict Method for Objects of Class modellist
#' 
#' Forecasting Bayesian VAR objects in a list of class 'modellist'.
#' 
#' @param object an object of class 'modellist'.
#' @param ... additional arguments.
#' 
#' @return A time-series object of class 'bvarprdlist'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 1:2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
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
#' model <- add_posterior_coefficients(model)
#' 
#' # Add forecast input data
#' model <- add_forecast_input(model, n_ahead = 10)
#' model <- add_posterior_forecasts(model)
#' 
#' # Generate forecasts
#' bvar_pred <- predict(model, n_ahead = 10)
#' 
#' 
#' @export
#' @method predict modellist
predict.modellist <- function(object, ...) {
  result <- lapply(object, predict, ...)
  class(result) <- c("bvarprdlist", "list")
  return(result)
}
