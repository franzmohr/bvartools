#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list of class 'expandingwindow', usually, the output of a call to
#' \code{\link{use_expanding_window}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'expandingwindow'.
#' 
#' @examples
#' 
#' data("us_macrodata")
#' 
#' # Starting period of the forecasting exercise
#' start_period <- 2007
#' 
#' # AR(1) models as benchmark
#' model <- create_bvarmodel(data = us_macrodata, p = 1,
#'                           deterministic = "none",
#'                           error = "gamma",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Obtain objects for expanding window estimation
#' model <- use_expanding_window(model, start = start_period)
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(shape = 3, rate = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @export
add_initial_values.expandingwindow <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- add_initial_values(object[[i]], ...)
  }
  
  return(object)
}