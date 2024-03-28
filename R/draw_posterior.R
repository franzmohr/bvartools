#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions. This is a generic function.
#' 
#' @param object a list of model specifications. Usually, the output of a call to 
#' \code{\link{create_var_model}} or \code{\link{create_vec_model}} in
#' combination with \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param ... arguments passed forward to method.
#' 
#' @export
draw_posterior <- function(object, ...){
  UseMethod("draw_posterior")
}
