#' Posterior Simulation for BVAR Models
#' 
#' Produces draws from the posterior distributions of Bayesian VAR models.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a
#' call to \code{\link{create_var_model}} in combination with
#' \code{\link{add_priors}} and \code{\link{add_initial_values}}.
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
#' The function also supports structural BVAR models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VAR model, the structural coefficients are drawn jointly with the other coefficients.
#' 
#' @return An object of class \code{"bvar"}.
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
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_var_model(e1, p = 2, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' object <- bvarpost(model)
#' 
#' @export
bvarpost <- function(object) {
  
  # Check if the input is suitable for the posterior simulation functions
  check_bvarpost_input(object)
  
  if (object[["model"]][["tvp"]]) {
    object <- .bvartvpalg(object)
  } else {
    object <- .bvaralg(object)
  }
  
  if (!is.null(object[["data"]][["z"]])) {
    if (object[["model"]][["tvp"]]) {
      for (i in 1:length(object[["posteriors"]][["a"]][["coeffs"]])) {
        object[["posteriors"]][["a"]][["coeffs"]][[i]]  <- coda::as.mcmc(object[["posteriors"]][["a"]][["coeffs"]][[i]])
      }
    } else {
      object[["posteriors"]][["a"]][["coeffs"]] <- coda::as.mcmc(object[["posteriors"]][["a"]][["coeffs"]])
    }
  }
  
  if (object[["model"]][["error"]] %in% c("sv", "sv-covar")) {
    for (i in 1:length(object[["posteriors"]][["sigma"]][["coeffs"]])) {
      object[["posteriors"]][["sigma"]][["coeffs"]][[i]]  <- coda::as.mcmc(object[["posteriors"]][["sigma"]][["coeffs"]][[i]])
    }
  } else {
    object[["posteriors"]][["sigma"]][["coeffs"]] <- coda::as.mcmc(object[["posteriors"]][["sigma"]][["coeffs"]])
  }
  
  
  return(object)
}