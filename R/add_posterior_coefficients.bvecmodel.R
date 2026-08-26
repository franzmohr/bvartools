#' Posterior Simulation of Model Coefficients
#' 
#' Forwards model input to posterior simulation functions for vector error correction models.
#' 
#' @param object an object of class 'bvecmodel', usually, a result of a
#' call to \code{\link{create_bvecmodel}} in combination with
#' \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to the model in argument \code{object}.
#' If \code{NULL} (default), internal functions are used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The function implements commonly used posterior simulation algorithms for Bayesian VAR models with
#' both constant and time varying parameters (TVP) as well as stochastic volatility. It can produce posterior
#' draws for standard BVAR models with independent normal-Wishart priors, which can be augmented by stochastic
#' search variable selection (SSVS) as proposed by Geroge et al. (2008) or Bayesian variable selection (BVS)
#' as proposed in Korobilis (2013). Both SSVS or BVS can also be applied to the covariances of the error term.
#' 
#' The implementation follows the descriptions in Chan et al. (2019), George et al. (2008) and Korobilis (2013).
#' For all approaches the SUR form of a VAR model is used to obtain posterior draws. The algorithm is implemented
#' in C++ to reduce calculation time.
#' 
#' The function also supports structural BVEC models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VEC model, the structural coefficients are drawn jointly with the other coefficients.
#' 
#' @return An object of class 'bvecmodel'.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
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
#' @export
add_posterior_coefficients.bvecmodel <- function(object, posterior_function = NULL, ...){
  
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
        .check_bvecpost_input(object)
        
        algorithm <- object[["model"]][["algorithm"]]
        
        if (algorithm %in% c("VecNormalGamma", "VecNormalWishart", "VecNormalStochvol",
                             "VecTvpGamma", "VecTvpWishart", "VecTvpStochvol")) {
          object <- switch(algorithm,
                           VecNormalGamma = .VecNormalGammaCoefficients(object),
                           VecNormalStochvol = .VecNormalStochvolCoefficients(object),
                           VecNormalWishart = .VecNormalWishartCoefficients(object),
                           VecTvpGamma = .VecTvpGammaCoefficients(object),
                           VecTvpStochvol = .VecTvpStochvolCoefficients(object),
                           VecTvpWishart = .VecTvpWishartCoefficients(object))
        } else {
          stop("Algorithm '", algorithm, "' not supported.")
        }
        
        for (i in c("a", "beta", "psi", "u_sigma_inv", "u_omega_inv")) {
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