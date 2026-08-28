#' Add Log-Likelihood
#' 
#' Forwards its elements to methods used to calculate posterior log-likelihoods.
#' 
#' @param object an object of class 'expandingwindow', usually, the result of a call to
#' \code{\link{add_posterior_coefficients}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'expandingwindow'.
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
#'                           iterations = 10, burnin = 10)
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
#' # Add log-likelihoods
#' model <- add_posterior_loglik(model)
#' 
#' @export
add_posterior_loglik.expandingwindow <- function(object, ...){
  
  orig_class <- class(object)
  object <- lapply(object, add_posterior_loglik, ...)
  class(object) <- orig_class
  
  return(object)
}