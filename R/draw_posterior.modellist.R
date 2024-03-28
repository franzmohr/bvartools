#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions.
#' 
#' @param object a list of models, usually, the output of a call to \code{\link{create_var_model}}
#' or \code{\link{create_vec_model}} in combination with \code{\link{add_priors}} and
#' \code{\link{add_initial_values}}.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- create_var_model(e1, p = 1:2, deterministic = 2,
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior.modellist <- function(object, mc.cores = NULL, ...){
  
  cat("Simulating...\n")
  
  if (is.null(mc.cores)) {
    object <- lapply(object, draw_posterior, ...)
  } else {
    object <- parallel::mclapply(object, draw_posterior, mc.cores = mc.cores, mc.preschedule = FALSE, ...)
  }
  
  class(object) <- append("posteriorlist", class(object))
  
  return(object)
}