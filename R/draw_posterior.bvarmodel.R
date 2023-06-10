#' Posterior Simulation
#' 
#' Forwards model input to posterior simulation functions.
#' 
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{gen_var}} in
#' combination with \code{\link{add_priors}}.
#' @param FUN the function to be applied to each model in argument \code{object}.
#' If \code{NULL} (default), the internal functions \code{\link{bvarpost}} is used.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return For multiple models a list of objects of class \code{bvarlist}.
#' For a single model the object has the class of the output of the applied posterior
#' simulation function. In case the package's own functions are used, this will
#' result in an object of class \code{"bvar"}.
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
draw_posterior.bvarmodel <- function(object, FUN = NULL, mc.cores = NULL, ...){
  
  # Check if it's only one model
  only_one_model <- FALSE
  if ("data" %in% names(object)) {
    object <- list(object)
    only_one_model <- TRUE
  }
  
  if (length(object) > 1) {
    cat("Estimating models...\n")
  } else {
    cat("Estimating model...\n") 
  }
  
  if (is.null(mc.cores)) {
    object <- lapply(object, .posterior_bvarmodel, use = FUN)
  } else {
    object <- parallel::mclapply(object, .posterior_bvarmodel, use = FUN,
                                 mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  
  if (only_one_model) {
    object <- object[[1]]
  } else {
    class(object) <- append("bvarlist", class(object)) 
  }
  
  return(object)
}

# Helper function to implement try() functionality
.posterior_bvarmodel <- function(object, use) {
  
  if (is.null(use)) {
    object <- try(bvarpost(object))
  } else {
    # Apply own function
    object <- try(use(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(object, list(error = TRUE))
  }
  
  return(object)
}
