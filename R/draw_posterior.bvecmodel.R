#' Posterior Simulation for Vector Error Correction Models
#' 
#' Forwards model input to posterior simulation functions for vector error correction models.
#' 
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{gen_vec}}
#' in combination with \code{\link{add_priors}}.
#' @param FUN the function to be applied to each list element in argument \code{object}.
#' If \code{NULL} (default), the internal function \code{\link{bvecpost}} is used.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return For multiple models a list of objects of class \code{bvarlist}.
#' For a single model the object has the class of the output of the applied posterior
#' simulation function. In case the package's own functions are used, this will
#' be \code{"bvec"}.
#' 
#' @references
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- gen_vec(e6, p = 1, r = 1, const = "restricted",
#'                  iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior.bvecmodel <- function(object, FUN = NULL, mc.cores = NULL, ...){
  
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
    object <- lapply(object, .posterior_bvecmodel, use = FUN)
  } else {
    object <- parallel::mclapply(object, .posterior_bvecmodel, use = FUN,
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
.posterior_bvecmodel <- function(object, use) {
  
  if (is.null(use)) {
    object <- try(bvecpost(object))
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
