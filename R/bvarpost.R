#' Posterior Simulation for BVAR Models
#' 
#' Produces draws from the posterior distributions of Bayesian VAR models.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' in combination with \code{\link{add_priors}}.
#' 
#' @details The function implements commonly used posterior simulation algorithms for Bayesian VAR models.
#' It can produce posterior draws for standard BVAR models with independent normal-Wishart priors, which can
#' be augmented by stochastic search variable selection (SSVS) as proposed by Geroge et al. (2008) or Bayesian
#' variable selection (BVS) as proposed in Korobilis (2013). Both SSVS and BVS can also be applied to the
#' covariances of the error term.
#' 
#' The implementation follows the description in Chan et al. (2019), George et al. (2008) and Korobilis (2013).
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
#' \url{https://doi.org/10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \url{https://doi.org/10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Get data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- gen_var(e1, p = 2, deterministic = "const",
#'                  iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#'
#' # Obtain posterior draws 
#' object <- bvarpost(model)
#' 
#' @export
bvarpost <- function(object) {
  
  object <- .bvaralg(object) # Use C++ code to draw posteriors
  
  A0 <- NULL
  if (object$model$structural) {
    k <- NCOL(object$data$Y)
    pos <- which(lower.tri(diag(1, k)))
    if ("list" %in% class(object$posteriors[["a0"]])) {
      if ("lambda" %in% names(object$posteriors[["a0"]])) {
        A0 <- list("coeffs" = matrix(diag(1, k), k * k, object$model$iterations),
                   "lambda" = matrix(diag(1, k), k * k, object$model$iterations))
        A0[["coeffs"]][pos, ] <- object$posteriors[["a0"]][["coeffs"]]
        A0[["lambda"]][pos, ] <- object$posteriors[["a0"]][["lambda"]]
      }
    } else {
      A0 <- matrix(diag(1, k), k * k, object$model$iterations)
      A0[pos, ] <- object$posteriors$a0 
    }
  }
  
  # Create bvar object
  object <- bvar(data = NULL, exogen = NULL,
                 y = object$data$Y, x = object$data$Z,
                 A0 = A0,
                 A = object$posteriors$a,
                 B = object$posteriors$b,
                 C = object$posteriors$c,
                 Sigma = object$posteriors$sigma)
  
  return(object)
}