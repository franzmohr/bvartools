#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a VAR model, which was produced by
#' function \code{\link{create_var_model}} in combination with \code{\link{add_priors}}.
#'
#' @param object list of class 'bvarmodel'.
#' @param method character specifying the method of how initial values are generated.
#' Default is \code{"ols"}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details 
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"ols"}}{Inital values are equal to estimates from an unrestricted LS regression.}
#'   \item{\code{"prior"}}{Initial values are drawn from the prior. Not possible for uninformative priors.}
#' }
#' 
#' In case \code{method = "ols"}, the initial draw of \eqn{a} is the result of
#' the unrestricted LS regression \eqn{(Z^{\prime}Z)^{-1}Z^{\prime}y} of the model
#' \eqn{y_{t} = Z_{t} a + u_{t}}, where \eqn{Z} is the matrix of regressors in
#' SUR form and \eqn{y} is the vector of endogenous variables.
#' 
#' The initial draw of \eqn{\Sigma^u} is the sum of squared
#' residuals divided by the number of observations, i.e. \eqn{\frac{uu^{\prime}}{T}},
#' with \eqn{u} as the residuals of the LS regression.
#' 
#' In case of a model with time varying parameters (TVP), the initial states are
#' obtained using the approach specified in argument \code{method}. However, the
#' initial draws of the error variances of the state equations are always drawn
#' from their prior distributions.
#' 
#' @return An object of class 'bvarmodel'.
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
#' @export
add_initial_values.bvarmodel <- function(object, method = "ols", ...){
  
  if (!"priors" %in% names(object)) {
    stop("No information on priors found in argument 'object'.")
  }
  
  if (!method %in% c("ols", "prior")) {
    stop("Argument 'method' can be 'ols' or 'prior' for BVAR models.")
  }
  
  y <- t(object[["data"]][["y"]])
  z <- object[["data"]][["z"]]
  k <- nrow(y)
  tt <- ncol(y)
  
  # Initial values from LS estimates ----
  if (method == "ols") {
    
    u <- y
    
    # Coefficients
    if (!is.null(z)) {
      
      p <- object[["model"]][["p"]]
      m <- object[["model"]][["m"]]
      s <- object[["model"]][["s"]]
      n <- object[["model"]][["n"]]
      nparams <- k * p + m * (s + 1) + n
      
      if (tt >= nparams) {
        ols <- solve(crossprod(z)) %*% crossprod(z, matrix(y))
        
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["a"]] <- matrix(ols, length(ols) * tt)
          object[["initial"]][["a_init"]] <- ols
        } else {
          object[["initial"]][["a"]] <- ols 
        }
        u <- matrix(matrix(y) - z %*% ols, NROW(y))
      } else {
        warning("Not enough observations for LS-based initial values. Setting initial values of coefficients to 0.")
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["a"]] <- matrix(0, ncol(z) * tt)
          object[["initial"]][["a_init"]] <- matrix(0, ncol(z))
        } else {
          object[["initial"]][["a"]] <- matrix(0, ncol(z))
        }
      }
    }
    
    # Covariances
    if (object$model$error %in% c("gamma-covar", "sv-covar") & k > 1) {
      y_covar <- kronecker(-t(u), diag(1, k))
      pos <- NULL
      for (j in 1:k) {pos <- c(pos, (j - 1) * k + 1:j)}
      y_covar <- y_covar[, -pos]
      psi <- solve(crossprod(y_covar)) %*% crossprod(y_covar, matrix(u))
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["psi"]] <- matrix(psi, length(psi) * tt)
        object[["initial"]][["psi_init"]] <- matrix(psi, length(psi))
      } else {
        object[["initial"]][["psi"]] <- psi
      }
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
      a <- a_mu + chol(a_vinv) %*% stats::rnorm(length(a_mu)) 
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["a"]] <- matrix(a, length(a) * tt)
        object[["initial"]][["a_init"]] <- matrix(a, length(a))
      } else {
        object[["initial"]][["a"]] <- a
      }
    }
    
    # Covariances
    if (!is.null(object$priors$psi) & k > 1) {
      psi_mu <- object$priors$psi$mu
      psi_vinv <- object$priors$psi$v_i
      psi <- psi_mu + chol(psi_vinv) %*% stats::rnorm(k * (k - 1) / 2)
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["psi"]] <- matrix(psi, length(psi) * tt)
        object[["initial"]][["psi_init"]] <- matrix(psi, length(psi))
      } else {
        object[["initial"]][["psi"]] <- psi
      }
    }
    
    u <- NULL
  }
  
  # Variances of state equations ----
  object <- .add_initial_values_state_errors(object)
  
  # Initial values for errors
  object <- .add_initial_values_measurement_errors(object = object,
                                                   method = method,
                                                   u = u)
  
  
  return(object)
}