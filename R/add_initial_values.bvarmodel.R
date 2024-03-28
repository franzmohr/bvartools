#' Add Initial Values to a Bayesian VAR Model
#'
#' Adds initial values to a VAR model, which was produced by
#' function \code{\link{create_var_model}} in combination with \code{\link{add_priors}}.
#'
#' @param object a list, usually, the output of a call to \code{\link{create_var_model}}.
#' @param method a character specifying the method of how initial values are generated.
#' Defaults to \code{"ols"}, which is not entirely available for TVP and SV models. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details 
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"ols"}}{Inital values are equal to estimates from unrestricted LS regressions.}
#'   \item{\code{"prior"}}{Initial values are drawn from the prior. Not possible for uninformative priors.}
#' }
#' 
#' In case \code{method = "ols"}, the initial draw of \eqn{a} is the result of
#' the unrestricted LS regression \eqn{(Z^{\prime}Z)^{-1}Z^{\prime}y}, where \eqn{Z}
#' is the matrix of regressors in SUR form and \eqn{y} is the vector of endogenous
#' variables.
#' The initial draw of \eqn{V} is the sum of squared
#' residuals divided by the number of observations, i.e. \eqn{\frac{vv^{\prime}}{T}},
#' with \eqn{v} as the residuals of the LS regression.
#' 
#' @return An object of class 'bvarmodel'.
#' 
#' @examples 
#' 
#' # Prepare data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- create_var_model(e1, p = 2, deterministic = 2,
#'                           iterations = 100, burnin = 10)
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @export
add_initial_values.bvarmodel <- function(object, method = "ols", ...){
  
  if (!"priors" %in% names(object)) {
    stop("No information on priors found in argument 'object'.")
  }
  
  y <- t(object[["data"]][["Y"]])
  z <- object[["data"]][["SUR"]]
  k <- nrow(y)
  tt <- ncol(y)
  if (!is.null(z)) {
    est_vars <- ncol(z)
  } else {
    est_vars <- 0
  }
  
  # Initial values from LS estimates ----
  if (method == "ols") {
    
    u <- y
    
    # Coefficients
    if (!is.null(z)) {
      
      if (nrow(z) / tt == k) {
        ols <- solve(crossprod(z)) %*% crossprod(z, matrix(y))
        object[["initial"]][["a"]] <- ols
        u <- matrix(matrix(y) - z %*% ols, NROW(y))
      } else {
        warning("Not enough observations for LS-based initial values. Setting initial values of coefficients to 0.")
        object[["initial"]][["a"]] <- matrix(0, ncol(z))
      }
      
      # Time varying parameters
      if (object[["model"]][["tvp"]]) {
        stop("Initial values for TVP not implemented yet.")
        object[["initial"]][["coefficients"]][["sigma_i"]] <- diag(c(1 / object[["priors"]][["coefficients"]][["rate"]]), ncol(z))
      }
    }
    
    # Covariances
    if (!is.null(object$priors$psi)) {
      y_covar <- kronecker(-t(u), diag(1, k))
      pos <- NULL
      for (j in 1:k) {pos <- c(pos, (j - 1) * k + 1:j)}
      y_covar <- y_covar[, -pos]
      psi <- solve(crossprod(y_covar)) %*% crossprod(y_covar, matrix(u))
      object[["initial"]][["psi"]] <- psi
      Psi <- diag(1, k)
      for (j in 2:k) {
        Psi[j, 1:(j - 1)] <- t(psi[((j - 2) * (j - 1) / 2) + 1:(j - 1), 1])
      }
      u <- Psi %*% u
    }
  }
  
  # Initial values from priors ----
  if (method == "prior") {
    # Coefficients
    if (!is.null(z)) {
      a_mu <- object$priors$a$mu
      a_vinv <- object$priors$a$v_i
      object[["initial"]][["a"]] <- a_mu + chol(a_vinv) %*% stats::rnorm(ncol(z))
      
      # Time varying parameters
      if (object[["model"]][["tvp"]] & !is.null(z)) {
        stop("Initial values for TVP not implemented yet.")
        object[["initial"]][["coefficients"]][["sigma_i"]] <- diag(c(1 / object[["priors"]][["coefficients"]][["rate"]]), ncol(z))
      }
    }
    
    # Covariances
    if (!is.null(object$priors$psi) & k > 1) {
      psi_mu <- object$priors$psi$mu
      psi_vinv <- object$priors$psi$v_i
      object[["initial"]][["psi"]] <- psi_mu + chol(psi_vinv) %*% stats::rnorm(k * (k - 1) / 2)
    }
    
    u <- NULL
  }
  
  # Initial values for errors
  object <- .add_initial_values_measurement_errors(object = object, method = method, u = u)
  
  
  return(object)
}