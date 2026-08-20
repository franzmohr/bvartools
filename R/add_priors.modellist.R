#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list of class 'modellist', usually, the output of a call
#' to \code{\link{create_bvarmodel}} or \code{\link{create_bvecmodel}}.
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
#' # Create model
#' model <- create_bvarmodel(e1, p = 0:2, deterministic = "const",
#'                           iterations = 10, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' @export
add_priors.modellist <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- add_priors(object[[i]], ...)
  }
  
  return(object)
}