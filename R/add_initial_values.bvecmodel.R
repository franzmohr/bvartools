#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a VEC model, which was produced by
#' function \code{\link{create_vec_model}} in combination with \code{\link{add_priors}}.
#'
#' @param object list of class 'bvecmodel'.
#' @param method character specifying the method of how initial values are generated.
#' Defaults is \code{"maxlik"}. Different approaches are used for TVP and SV models. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"maxlik"}}{Inital values are equal to estimates from maximum likelihood regressions.}
#'   \item{\code{"prior"}}{Initial values are drawn from the prior. Not possible for uninformative priors.}
#' }
#' 
#' In case \code{method = "maxlik"}, the initial draw of \eqn{a} is the result of
#' a maximum likelihood estimation of the reduced rank model. When used with a
#' Wishart prior, the initial draw of \eqn{\Sigma^u} is the sum of squared
#' residuals of ML regression divided by the number of observations, i.e.
#' \eqn{\frac{uu^{\prime}}{T}}. In all other cases, the diagonal elements of
#' \eqn{\Sigma^u} are set to the variances of the variables in \eqn{u}.
#' 
#' In case \code{method = "prior"}, all initial draws in the model are random
#' draws from the respective prior distributions.
#' 
#' In case of a model with time varying parameters (TVP), the initial states are
#' obtained using the approach specified in argument \code{method}. However, the
#' initial draws of the error variances of the state equations are always drawn
#' from their prior distributions. The value of the autocorrelation coefficients
#' \code{rho} of that states of \eqn{\beta} are set to the value in the prior
#' specification.
#' 
#' @return An object of class 'bvecmodel'.
#' 
#' @examples
#' 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_vec_model(e6, p = 4, r = 1,
#'                           const = "unrestricted", seasonal = "unrestricted",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @export
add_initial_values.bvecmodel <- function(object, method = "maxlik", ...){
  
  if (!"priors" %in% names(object)) {
    stop("No information on priors found in argument 'object'.")
  }
  
  if (!method %in% c("maxlik", "prior")) {
    stop("Invalid specification of argument 'method'.")
    stop("Argument 'method' can be 'maxlik' or 'prior' for BVEC models.")
  }
  
  y <- t(object[["data"]][["y"]])
  w <- t(object[["data"]][["w"]])
  if (!is.null(object[["data"]][["x"]])) {
    x <- t(object[["data"]][["x"]]) 
  } else {
    x <- NULL
  }
  z <- object[["data"]][["z"]]
  r <- object[["model"]][["rank"]]
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n_unrestricted"]]
  
  k_ect <- ncol(object[["data"]][["w"]])
  n_alpha <- k * r
  n_beta <- k_ect * r
  n_ect <- k_ect * k
  tt <- ncol(y)
  if (is.null(z)) {
    n_z <- 0
  } else {
    n_z <- ncol(z)
  }
  
  # Maximum likelihood ----
  if (method == "maxlik") {
    
    u <- y
    
    ## Coefficients ----
    if (tt >= r + k * (p - 1) + m * s + n) {
      
      if (r > 0) {
        
        if (k * (p - 1) + m * s + n > 0) {
          M <- diag(tt) - crossprod(x, solve(tcrossprod(x))) %*% x
          R0 <- y %*% M # Residuals of regression of y on x
          R1 <- w %*% M # Residuals of regression of w on x 
        } else {
          R0 <- y
          R1 <- w
        }
        S00_inv <- solve(tcrossprod(R0) / tt)
        S01 <- tcrossprod(R0, R1) / tt
        S10 <- tcrossprod(R1, R0) / tt
        S11 <- tcrossprod(R1) / tt
        S11_sqrt_inv <- solve(.mroot(S11))
        lambda <- eigen(S11_sqrt_inv %*% S10 %*% S00_inv %*% S01 %*% t(S11_sqrt_inv))#, symmetric = TRUE)
        
        beta <- t(crossprod(matrix(lambda$vectors[, 1:r] , nrow(w)), S11_sqrt_inv))
        z[, 1:n_alpha] <- kronecker(t(crossprod(beta, w)), diag(1, k))
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["beta"]] <- matrix(beta, length(beta) * tt)
          object[["initial"]][["beta_init"]] <- matrix(beta)
        } else {
          object[["initial"]][["beta"]] <- beta 
        }
      }
      
      if (n_z > 0) {
        ml <- solve(crossprod(z)) %*% crossprod(z, matrix(y))
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["a"]] <- matrix(ml, length(ml) * tt)
          object[["initial"]][["a_init"]] <- ml
        } else {
          object[["initial"]][["a"]] <- ml 
        }
        u <- matrix(matrix(y) - z %*% ml, NROW(y)) 
      }
      
    } else {
      warning("Not enough observations for ML-based initial values. Setting initial values of coefficients to 0.")
      if (n_z > 0) {
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["a"]] <- matrix(0, ncol(z) * tt)
          object[["initial"]][["a_init"]] <- matrix(0, ncol(z))
        } else {
          object[["initial"]][["a"]] <- matrix(0, ncol(z))
        } 
      }
    }
    
    ## Covariances ----
    if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & k > 1) {
      y_covar <- kronecker(-t(u), diag(1, k))
      pos <- NULL
      for (j in 1:k) {pos <- c(pos, (j - 1) * k + 1:j)}
      y_covar <- y_covar[, -pos]
      psi <- solve(crossprod(y_covar)) %*% crossprod(y_covar, matrix(u))
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["psi"]] <- matrix(psi, length(psi) * tt)
        object[["initial"]][["psi_init"]] <- psi
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
    if (!is.null(object$priors$a)) {
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
    
    if (r > 0) {
      beta <- matrix(0, n_ect / k, object$model$rank)
      beta[1:object$model$rank, 1:object$model$rank] <- diag(1, object$model$rank)
      object[["initial"]][["beta"]] <- beta
      z[, 1:n_alpha] <- kronecker(t(crossprod(beta, w)), diag(1, k))
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["beta"]] <- matrix(beta, length(beta) * tt)
        object[["initial"]][["beta_init"]] <- matrix(beta)
      } else {
        object[["initial"]][["beta"]] <- beta
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
  
  if (n_z > 0) {
    object[["data"]][["z"]] <- z 
  }
  
  # Variances of state equations ----
 object <- .add_initial_values_state_errors(object)
  
  # Initial values for errors
  object <- .add_initial_values_measurement_errors(object = object,
                                                   method = method,
                                                   u = u)
  
  return(object)
}

# Square root of a matrix
.mroot <- function(M){
  eig <- eigen(M)
  if (length(eig$values) == 1){
    val <- matrix(sqrt(eig$values), 1)
  } else {
    val <- diag(sqrt(eig$values))
  }
  R <- eig$vectors %*% val %*% t(eig$vectors)
  return(R)
}
