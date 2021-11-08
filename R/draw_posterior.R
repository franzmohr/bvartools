#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions. This is a generic function.
#' 
#' @param object a list of model specifications. Usually, the output of a call to 
#' \code{\link{gen_var}}, \code{\link{gen_vec}} or \code{\link{gen_dfm}} in combination with \code{\link{add_priors}}.
#' @param ... arguments passed forward to method.
#' 
#' @export
draw_posterior <- function(object, ...){
  UseMethod("draw_posterior")
}
