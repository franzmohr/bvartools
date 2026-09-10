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
  .transition_message("draw_posterior", "add_posterior_coefficients",
                      note = paste("Forecasts and the log likelihood are added separately, by",
                                   "'add_posterior_forecasts()' and 'add_posterior_loglik()'."),
                      also = c("bvarpost", "bvecpost", "dfmpost"))

  UseMethod("draw_posterior")
}
