#' Add Log-Likelihood
#' 
#' Calculates and adds posterior log-likelihoods to an object of class 'bvecmodel'.
#' 
#' @param object an object of class 'bvecmodel', usually, the result of a call to
#' \code{\link{add_posterior_coefficients}}.
#' @param ... additional arguments.
#' 
#' @return An object of class 'bvecmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_bvecmodel(e6, p = 4, r = 1,
#'                           const = "unrestricted",
#'                           seasonal = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Add log-likelihoods
#' model <- add_posterior_loglik(model)
#' 
#' 
#' @export
#' @method add_posterior_loglik bvecmodel
add_posterior_loglik.bvecmodel <- function(object, ...) {
  
  # Input checks
  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Object does not contain posterior draws in posteriors$u_sigma_inv.")
  }
  
  algorithm <- object[["model"]][["algorithm"]]
  
  if (algorithm %in% c("VecNormalWishart")) {
    object <- switch(algorithm,
                     VecNormalGamma = .VecNormalGammaLogLik(object),
                     VecNormalStochvol = .VecNormalStochvolLogLik(object),
                     VecNormalWishart = .VecNormalWishartLogLik(object),
                     VecTvpGamma = .VecTvpGammaLogLik(object),
                     VecTvpStochvol = .VecTvpStochvolLogLik(object),
                     VecTvpWishart = .VecTvpWishartLogLik(object))
  } else {
    stop("Algorithm not implemented yet.")
  }
  
  # Calculate log likelihoods
  temp_mcpar <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["loglik"]] <- coda::mcmc(object[["posterior"]][["loglik"]], start = temp_mcpar[1], end = temp_mcpar[2], thin = temp_mcpar[3])
  
  return(object)
}
