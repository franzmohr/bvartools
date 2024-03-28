#' Posterior Simulation for Vector Autoregressive Models
#' 
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#' 
#' @param object a list of models, usually, the output of a call to
#' \code{\link{create_var_model}} in combination with \code{\link{add_priors}}
#' and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to each model in argument \code{object}.
#' If \code{NULL} (default), the internal functions \code{\link{bvarpost}} is used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of the class of the output of the applied posterior
#' simulation function. In case the package's own function is used, this will
#' result in an object of class \code{"bvar"}.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- create_var_model(e1, p = 2, deterministic = 2,
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
draw_posterior.bvarmodel <- function(object, posterior_function = NULL, ...){

  if (is.null(posterior_function)) {
    object <- try(bvarpost(object))
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(object, list(error = TRUE))
  }
  
  return(object)
}