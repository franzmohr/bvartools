#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list, usually, the output of a call to
#' \code{\link{create_var_model}} or \code{\link{create_vec_model}} in
#' combination with \code{\link{add_priors}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of models.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_var_model(e1, p = 0:2, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
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