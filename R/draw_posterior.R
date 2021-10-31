#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions.
#' 
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{gen_var}},
#' \code{\link{gen_vec}} or \code{\link{gen_dfm}} in combination with \code{\link{add_priors}}.
#' @param ... arguments passed forward to method.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- gen_var(e1, p = 1:2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior <- function(object, ...){
  UseMethod("draw_posterior")
}
