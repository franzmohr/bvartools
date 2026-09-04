#' Stochastic Search Variable Selection Prior
#' 
#' Calculates the priors for a Bayesian VAR model, which employs stochastic search variable selection (SSVS).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{create_bvarmodel}}
#' or \code{\link{create_bvecmodel}}.
#' @param tau a numeric vector of two elements containing the prior standard errors of restricted
#' variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#' as its second. Default is \code{c(0.05, 10)}.
#' @param semiautomatic an optional numeric vector of two elements containing the factors by which
#' the standard errors associated with an unconstrained least squares estimate of the VAR model are
#' multiplied to obtain the prior standard errors of restricted (\eqn{\tau_0}) and unrestricted
#' (\eqn{\tau_1}) variables. This is the semiautomatic approach described in George et al. (2008).
#' 
#' @return A list containing the vectors of prior standard deviations for restricted
#' and unrestricted variables, respectively.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' @examples
#' 
#' # Prepare data
#' data("e1")
#' train <- diff(log(e1))
#' 
#' # Generate model input
#' model <- create_bvarmodel(train, p = 2, deterministic = "const",
#'                           iterations = 5000, burnin = 1000)
#' 
#' # Obtain SSVS prior
#' prior <- ssvs_prior(model, semiautomatic = c(.1, 10))
#' 
#' @export
ssvs_prior.bvarmodel <- function(object, tau = c(0.05, 10), semiautomatic = NULL) {
  
  if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
    stop("SSVS cannot be used with models with stochastic volatility.")
  }
  
  if (object[["model"]][["tvp"]]) {
    stop("SSVS cannot be used with models with time varying parameter.")
  }
  
  if (!is.null(semiautomatic)) {
    if (!"numeric" %in% class(semiautomatic)) {
      stop("Argument 'semiautomatic' must be a numeric vector of length 2.")
    }
    if (length(semiautomatic) != 2) {
      stop("Argument 'semiautomatic' must be a numeric vector of length 2.")
    }
  }
  
  z <- object[["data"]][["train"]][["z"]]
  
  if (!is.null(z)) {
    
    if (!is.null(semiautomatic)) {
      
      # Semiautomatic approach
      if (!is.null(z)) {
        
        y <- matrix(t(object[["data"]][["train"]][["y"]]))
        k <- object[["model"]][["k"]]
        tt <- nrow(object[["data"]][["train"]][["y"]])
        
        ols <- solve(crossprod(z)) %*% crossprod(z, y)
        u <- matrix(y - z %*% ols, k)
        sigma_ols <- tcrossprod(u) / (tt - ncol(z) / k) # OLS error covariance matrix
        cov_ols <- solve(crossprod(z, kronecker(diag(1, tt), solve(sigma_ols))) %*% z)
        se_ols <- matrix(sqrt(diag(cov_ols))) # OLS standard errors
        
        tau0 <- se_ols * semiautomatic[1] # Prior if excluded
        tau1 <- se_ols * semiautomatic[2] # Prior if included
      } else {
        tau0 <- NULL
        tau1 <- NULL
      }
      
    } else {
      tot_par <- NCOL(z)
      tau0 <- matrix(rep(tau[1], tot_par))
      tau1 <- matrix(rep(tau[2], tot_par))
    }
    
    result <- list("tau0" = tau0,
                   "tau1" = tau1)
    
  } else {
    result <- NULL
  }
  
  return(result)
}