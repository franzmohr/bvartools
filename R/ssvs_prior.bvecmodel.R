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
#' data("e6")
#' 
#' # Generate model input
#' object <- create_bvecmodel(e6, r = 1)
#' 
#' # Obtain SSVS prior
#' prior <- ssvs_prior(object, semiautomatic = c(.1, 10))
#' 
#' @export
ssvs_prior.bvecmodel <- function(object, tau = c(0.05, 10), semiautomatic = NULL) {
  
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
  
  y <- t(object[["data"]][["y"]])
  w <- t(object[["data"]][["w"]])
  tt <- NCOL(y)
  k <- NROW(y)
  
  # The output must have the same length as the number of columns in Z
  
  if (!is.null(object$data$z)) {
    
    tau0 <- NULL
    tau1 <- NULL
    
    # alpha coefficients ----
    if (object$model$rank > 0) {
      tau0 <- rbind(tau0, matrix(1, k * object$model$rank))
      tau1 <- rbind(tau1, matrix(1, k * object$model$rank))
    }
    
    # Non-alpha coefficients ----
    if (!is.null(object[["data"]][["x"]])) {
      
      x <- t(object[["data"]][["x"]])
      
      if (!is.null(semiautomatic)) {
        
        ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
        u <- y - ols %*% x
        sigma_ols <- tcrossprod(u) / (tt - nrow(x)) # OLS error covariance matrix
        cov_ols <- kronecker(solve(tcrossprod(x)), sigma_ols) # Sqrt of diagonal elements are the t-ratios
        se_ols <- sqrt(diag(cov_ols)) # OLS standard errors 
        
        tau0 <- append(tau0, se_ols * semiautomatic[1]) # Prior if excluded
        tau1 <- append(tau0, se_ols * semiautomatic[2]) # Prior if included
        
      } else {
        tau0 <- append(tau0, rep(tau[1], k * nrow(x)))
        tau1 <- append(tau1, rep(tau[2], k * nrow(x)))
      }
    }
    
    if (object[["model"]][["structural"]]) {
      message("Semiautomatic approach to SSVS not available for structural coefficients yet. Using values of argument 'tau' instead.")
      tau0 <- append(tau0, rep(tau[1], k * (k - 1) / 2))
      tau1 <- append(tau1, rep(tau[2], k * (k - 1) / 2))
    }
    
    result <- list("tau0" = matrix(tau0),
                   "tau1" = matrix(tau1))
    
  } else {
    result <- NULL
  }
  
  return(result)
}