#' Add Forecasts
#'
#' Calculates and adds forecasts to the elements of an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow', usually, the result of a call
#' to \code{\link{add_posterior_coefficients}} and \code{\link{add_forecast_input}}.
#' @param ... arguments passed forward to method.
#' @inheritParams add_posterior_coefficients.modellist
#' @inheritSection add_posterior_coefficients.modellist Parallel simulation
#' 
#' @return The object in \code{object} with forecast draws added to each of its models, as
#' described in \code{\link{add_posterior_forecasts.bvarmodel}}.
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
#'                           iterations = 10, burnin = 2)
#' # Number of iterations and burnin should be much higher.
#' 
#' model <- use_expanding_window(model, start = 1982.25)
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
#' @export
add_posterior_forecasts.expandingwindow <- function(object, ..., cores = 1){

  if (.use_cluster(object, cores)) {
    return(.simulate_models_in_parallel(object, add_posterior_forecasts, cores, ...))
  }
  
  orig_class <- class(object)
  object <- lapply(object, add_posterior_forecasts, ...)
  class(object) <- orig_class
  
  return(object)
}
