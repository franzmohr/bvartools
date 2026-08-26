#' Add Log-Likelihood
#' 
#' Calculates and adds posterior log-likelihoods to an object of class 'bvarmodel'.
#' 
#' @param object an object of class 'bvarmodel', usually, the result of a call to
#' \code{\link{add_posterior_coefficients}}.
#' @param ... additional arguments.
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
#' object <- add_posterior_coefficients(model)
#' 
#' # Add log-likelihoods
#' object <- add_posterior_loglik(object)
#' 
#' 
#' @export
#' @method add_posterior_loglik bvarmodel
add_posterior_loglik.bvarmodel <- function(object, ...) {
  
  # Input checks
  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Object does not contain posterior draws in posteriors$u_sigma_inv.")
  }
  
  algorithm <- object[["model"]][["algorithm"]]
  
  if (algorithm %in% c("VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                       "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")) {
    object <- switch(algorithm,
                     VarNormalGamma = .VarNormalGammaLogLik(object),
                     VarNormalStochvol = .VarNormalStochvolLogLik(object),
                     VarNormalWishart = .VarNormalWishartLogLik(object),
                     VarTvpGamma = .VarTvpGammaLogLik(object),
                     VarTvpStochvol = .VarTvpStochvolLogLik(object),
                     VarTvpWishart = .VarTvpWishartLogLik(object))
  } else {
    stop("Algorithm not implemented yet.")
  }
  
  # Calculate log likelihoods
  temp_mcpar <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["loglik"]] <- coda::mcmc(object[["posterior"]][["loglik"]], start = temp_mcpar[1], end = temp_mcpar[2], thin = temp_mcpar[3])
  
  return(object)
}
