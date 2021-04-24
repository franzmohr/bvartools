#' Posterior Simulation for Dynamic Factor Models
#' 
#' Produces draws from the posterior distributions of Bayesian dynamic factor models.
#' 
#' @param object an object of class \code{"dfmodel"}, usually, a result of a call to \code{\link{gen_dfm}}
#' in combination with \code{\link{add_priors}}.
#' 
#' @details The function implements the posterior simulation algorithm for Bayesian dynamic factor models.
#' 
#' The implementation follows the description in Chan et al. (2019) and C++ is used to reduce calculation time.
#' 
#' @return An object of class \code{"dfm"}.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' @examples
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(v_i = .01),
#'                     sigma_u = list(shape = 5, rate = 4),
#'                     a = list(v_i = .01),
#'                     sigma_v = list(shape = 5, rate = 4))
#' 
#' # Obtain posterior draws
#' object <- dfmpost(model)
#' 
#' @export
dfmpost <- function(object) {
  
  object <- .dfmalg(object) # Use C++ code to draw posteriors

  #x = object$data$X; lambda = object$posteriors$lambda; fac = object$posteriors$factor; sigma_v = object$posteriors$sigma_v; a = object$posteriors$a; sigma_u = object$posteriors$sigma_u
  
  object <- dfm(x = object$data$X,
                lambda = object$posteriors$lambda,
                fac = object$posteriors$factor,
                sigma_u = object$posteriors$sigma_u,
                a = object$posteriors$a,
                sigma_v = object$posteriors$sigma_v)
  
  return(object)
}