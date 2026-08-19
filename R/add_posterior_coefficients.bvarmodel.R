#' Posterior Simulation of Model Coefficients
#' 
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#' 
#' @param object an object of class 'bvarmodel', usually, a result of a
#' call to \code{\link{create_bvarmodel}} in combination with
#' \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to the model in argument \code{object}.
#' If \code{NULL} (default), internalal functions are used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details Unless \code{posterior_function} is specified, the function forwards
#' the model input to the package's own posterior functions.
#' 
#' @return An object of class 'bvarmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' @export
add_posterior_coefficients.bvarmodel <- function(object, posterior_function = NULL, ...){
  
  # This  allows to employ the method with other compatible classes
  class_of_object <- class(object)
  
  # Copy in case the simulation fails
  model <- object
  if ("posteriors" %in% names(model)) {
    model[["posteriors"]] <- NULL
  }
  
  if (is.null(posterior_function)) {
    object <- try(
      {
        # Check if the input is suitable for the posterior simulation functions
        .check_bvarpost_input(object)
        
        algorithm <- object[["model"]][["algorithm"]]
        
        if (algorithm %in% c("VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                             "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")) {
          object <- switch(algorithm,
                           VarNormalGamma = .VarNormalGammaCoefficients(object),
                           VarNormalStochvol = .VarNormalStochvolCoefficients(object),
                           VarNormalWishart = .VarNormalWishartCoefficients(object),
                           VarTvpGamma = .VarTvpGammaCoefficients(object),
                           VarTvpStochvol = .VarTvpStochvolCoefficients(object),
                           VarTvpWishart = .VarTvpWishartCoefficients(object))
        } else {
          stop("Algorithm '", algorithm, "' not supported.")
        }
        
        for (i in c("a", "psi", "u_sigma_inv", "u_omega_inv")) {
          if (!is.null(object[["posterior"]][[i]][["coeffs"]])) {
            object[["posterior"]][[i]][["coeffs"]] <- coda::as.mcmc(object[["posterior"]][[i]][["coeffs"]])
          }
          if (!is.null(object[["posterior"]][[i]][["lambda"]])) {
            object[["posterior"]][[i]][["lambda"]] <- coda::as.mcmc(object[["posterior"]][[i]][["lambda"]])
          }
          if (!is.null(object[["posterior"]][[i]][["sigma"]])) {
            object[["posterior"]][[i]][["sigma"]] <- coda::as.mcmc(object[["posterior"]][[i]][["sigma"]])
          }
        }
        
        object
      }
    )
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  }
  
  class(object) <- class_of_object
  
  return(object)
}