#' Stochastic Search Variable Selection Prior
#' 
#' Calculates the priors for a Bayesian VAR model, which employs stochastic search variable selection (SSVS).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' or \code{\link{gen_vec}}.
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
#' data <- diff(log(e1))
#' 
#' # Generate model input
#' object <- gen_var(data)
#' 
#' # Obtain SSVS prior
#' prior <- ssvs_prior(object, semiautomatic = c(.1, 10))
#' 
#' @export
ssvs_prior.bvarmodel <- function(object, tau = c(0.05, 10), semiautomatic = NULL) {
  
  if (object[["model"]][["error"]] %in% c("sv", "sv-covar")) {
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
  
  y <- t(object$data$y)
  tt <- NCOL(y)
  k <- NROW(y)
  
  if (!is.null(object$data$x)) {
    x <- t(object$data$x)
  } else {
    x <- NULL
  }
  z <- object$data$z
  
  if (!is.null(z)) {
    
    if (!is.null(semiautomatic)) {
      
      # Semiautomatic approach
      if (!is.null(x)) {
        ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
        u <- y - ols %*% x
        sigma_ols <- tcrossprod(u) / (tt - nrow(x)) # OLS error covariance matrix
        cov_ols <- kronecker(solve(tcrossprod(x)), sigma_ols) # Sqrt of diagonal elements are the t-ratios
        se_ols <- matrix(sqrt(diag(cov_ols))) # OLS standard errors 
        
        tau0 <- se_ols * semiautomatic[1] # Prior if excluded
        tau1 <- se_ols * semiautomatic[2] # Prior if included
      } else {
        tau0 <- NULL
        tau1 <- NULL
      }
      
      if (object[["model"]][["structural"]]) {
        message("Semiautomatic approach to SSVS not available for structural coefficients yet. Using values of argument 'tau' instead.")
        tau0 <- rbind(tau0, matrix(tau[1], k * (k - 1) / 2))
        tau1 <- rbind(tau1, matrix(tau[2], k * (k - 1) / 2))
      }
      
      if (object[["model"]][["error"]] %in% c("gamma-covar")) {
        message("Semiautomatic approach to SSVS not available for covariance coefficients yet. Using values of argument 'tau' instead.")
        tau0 <- rbind(tau0, matrix(tau[1], k * (k - 1) / 2))
        tau1 <- rbind(tau1, matrix(tau[2], k * (k - 1) / 2))
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