#' Add Initial Values to a List of Bayesian Models
#'
#' Adds initial values to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list, usually, the output of a call to
#' \code{\link{create_var_model}} or \code{\link{create_vec_model}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of models.
#' 
#' @examples 
#' 
#' # Prepare data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate models
#' model <- create_var_model(e1, p = 1:2, deterministic = 2,
#'                           iterations = 100, burnin = 10)
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @export
add_initial_values.modellist <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- add_initial_values(object[[i]], ...)
  }
  
  return(object)
}