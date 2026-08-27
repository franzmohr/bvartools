#' Model Selection Criteria
#'
#' Calculates model selection criteria for a list of Bayesian models.
#'
#' @param object an object of class 'modellist'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'selcritlist'.
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
#' model <- create_bvarmodel(data = train, p = 0:2, deterministic = "const",
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
#' # Calculate selection criteria
#' sel <- selection_criteria(model)
#' sel
#' 
#' @export
#' @method selection_criteria modellist
selection_criteria.modellist <- function(object, ...){
  
  object <- lapply(object, selection_criteria, ...)
  
  class(object) <- append("selcritlist", class(object))
  
  return(object)
}