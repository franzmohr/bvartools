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
ssvs_prior <- function(object, tau = c(0.05, 10), semiautomatic = NULL) {
  
  y <- t(object$data$Y)
  if (object$model$type == "VAR") {
    x <- t(object$data$Z)
  }
  if (object$model$type == "VEC") {
    if (!is.null(object$data$X)) {
      if (object$model$rank == 0) {
        x <- t(object$data$X)
      } else {
        x <- t(cbind(object$data$W, object$data$X)) 
      }
    } else {
      x <- t(object$data$W)
    }
  }
  
  tt <- NCOL(y)
  k <- NROW(y)
  tot_par <- k * NROW(x)
  if (object$model$structural) {
    tot_par <- tot_par + k * (k - 1) / 2
  }
  
  if (!is.null(semiautomatic)) {
    
    # Semiautomatic approach
    ols <- tcrossprod(y, x) %*% solve(tcrossprod(x)) # OLS estimates
    sigma_ols <- tcrossprod(y - ols %*% x) / (tt - nrow(x)) # OLS error covariance matrix
    cov_ols <- kronecker(solve(tcrossprod(x)), sigma_ols) # Sqrt of diagonal elements are the t-ratios
    se_ols <- matrix(sqrt(diag(cov_ols))) # OLS standard errors
    
    tau0 <- se_ols * semiautomatic[1] # Prior if excluded
    tau1 <- se_ols * semiautomatic[2] # Prior if included
    
    if (object$model$structural) {
      warning("Semiautomatic approach not available for structural variables. Using values of argument 'tau' instead.")
      
      tau0 <- rbind(tau0, matrix(tau[1], k * (k - 1) / 2))
      tau1 <- rbind(tau1, matrix(tau[2], k * (k - 1) / 2))
    }
    
    # Add structural variables here
  } else {
    tau0 <- matrix(rep(tau[1], tot_par))
    tau1 <- matrix(rep(tau[2], tot_par))
  }
  
  result <- list("tau0" = tau0,
                 "tau1" = tau1)
  
  return(result)
}