#' Posterior Simulation of Model Coefficients
#' 
#' Forwards model input to posterior simulation functions.
#' 
#' @param object an object of class 'modellist'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'modellist'.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- create_bvarmodel(e1, p = 1:2, deterministic = 2,
#'                           iterations = 20, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
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
#' @export
add_posterior_coefficients.modellist <- function(object, ...){
  
  object <- lapply(object, add_posterior_coefficients, ...)
  
  class(object) <- list("modellist", "list")
  
  return(object)
}