#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a VEC model, which was produced by
#' function \code{\link{create_bvecmodel}} in combination with \code{\link{add_priors}}.
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
#' obtained using the approach specified in argument \code{method}. The precisions
#' of the state equations start at the means of their gamma priors in case
#' \code{method = "maxlik"} and are drawn from those priors in case
#' \code{method = "prior"}. The autocorrelation coefficient \code{rho} of
#' the state equation of \eqn{\beta} is taken from the prior specification,
#' whether it is held there or drawn from the prior on it that
#' \code{coint$rho_min} and \code{coint$rho_max} set up.
#' 
#' @return The object in \code{object} with the element \code{initial} added, a list with
#' the starting values \code{beta} of the cointegration coefficients, a
#' \eqn{K_\beta r \times 1} matrix, or \eqn{T K_\beta r \times 1} with initial state
#' \code{beta_init} for time varying cointegration, and the starting values of the
#' remaining coefficients and the error term with the elements described in
#' \code{\link{add_initial_values.bvarmodel}}.
#'
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_bvecmodel(e6, p = 4, r = 1,
#'                           const = "unrestricted",
#'                           seasonal = "unrestricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' @family model set-up
#' @export
add_initial_values.bvecmodel <- function(object, method = "maxlik", ...){
  
  if (!"priors" %in% names(object)) {
    stop("No information on priors found in argument 'object'.")
  }
  
  if (!method %in% c("maxlik", "prior")) {
    stop("Argument 'method' can be 'maxlik' or 'prior' for BVEC models.")
  }
  
  y <- object[["data"]][["train"]][["y"]]
  w <- object[["data"]][["train"]][["w"]]
  if (!is.null(object[["data"]][["train"]][["x"]])) {
    x <- object[["data"]][["train"]][["x"]]
  } else {
    x <- NULL
  }
  z <- object[["data"]][["train"]][["z"]]
  r <- object[["model"]][["rank"]]
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  tt <- nrow(y)
  
  k_ect <- ncol(w)
  n_alpha <- k * r
  n_beta <- k_ect * r
  n_ect <- k_ect * k
  
  if (is.null(z)) {
    n_z <- 0
  } else {
    n_z <- ncol(z)
  }
  
  # 'y' arrives one row per period. The residual is carried as k x T, one column
  # per period, which is what the covariance block below and the error helpers
  # read it as, so it is transposed rather than reshaped: matrix(y, k) on the
  # wide series would fill each column with consecutive observations of the same
  # variable instead of one period across variables.
  u <- matrix(t(y), k)

  # Maximum likelihood ----
  if (method == "maxlik") {
    
    ## Coefficients ----
    if (tt >= r + k * (p - 1) + m * s + n) {

      # Outside the rank block below, because the least squares fit further down
      # stacks 'y' with matrix(y) whether or not there is a cointegration term
      # to estimate first. Transposing only under r > 0 meant a rank zero model
      # was fitted against the series stacked variable by variable rather than
      # period by period -- the wrong response for the SUR regressors, and wrong
      # in silence.
      y <- t(y)

      if (r > 0) {

        w <- t(w)

        beta <- .coint_ml(object)[["beta"]]
        ect <- crossprod(beta, w)
        z[, 1:n_alpha] <- kronecker(t(ect), diag(1, k))
        if (object[["model"]][["tvp"]]) {
          object[["initial"]][["beta"]] <- matrix(beta, length(beta) * tt)
          object[["initial"]][["beta_init"]] <- matrix(beta)
        } else {
          object[["initial"]][["beta"]] <- matrix(beta)
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
        u <- matrix(matrix(y) - z %*% ml, k) 
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

    # Both, as the maximum likelihood branch above does: 'y' arrives one row per
    # period, while the residual below is formed against the SUR regressors and
    # so needs the series stacked variable-within-period. Transposing only 'w'
    # left 'y' in its wide shape and the subtraction non-conformable.
    y <- t(y)
    w <- t(w)

    # Coefficients
    if (n_z > 0) {
      a_mu <- object[["priors"]][["a"]][["mu"]]
      a_vinv <- object[["priors"]][["a"]][["v_inv"]]
      if (all(diag(a_vinv) == 0)) {
        stop("All diagonal elements of the prior precision matrix of 'a' are zero.")
      }
      a <- a_mu + chol(a_vinv) %*% stats::rnorm(length(a_mu))
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["a"]] <- matrix(a, length(a_mu) * tt)
        object[["initial"]][["a_init"]] <- matrix(a, length(a_mu))
      } else {
        object[["initial"]][["a"]] <- a
      }
    }
    
    if (r > 0) {
      beta <- matrix(0, n_ect / k, object[["model"]][["rank"]])
      beta[1:object[["model"]][["rank"]], 1:object[["model"]][["rank"]]] <- diag(1, object[["model"]][["rank"]])
      object[["initial"]][["beta"]] <- beta
      z[, 1:n_alpha] <- kronecker(t(crossprod(beta, w)), diag(1, k))
      if (object[["model"]][["tvp"]]) {
        object[["initial"]][["beta"]] <- matrix(beta, length(beta) * tt)
        object[["initial"]][["beta_init"]] <- matrix(beta)
      } else {
        object[["initial"]][["beta"]] <- matrix(beta)
      }
      
    }
    
    if (n_z > 0) {
      u <- matrix(matrix(y) - z %*% a, k)
    }

    # Covariances
    if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & k > 1) {
      psi_mu <- object[["priors"]][["psi"]][["mu"]]
      psi_vinv <- object[["priors"]][["psi"]][["v_inv"]]
      psi <- psi_mu + chol(psi_vinv) %*% stats::rnorm(k * (k - 1) / 2)
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
  
  if (n_z > 0) {
    object[["data"]][["train"]][["z"]] <- z
  }

  # Inclusion indicators ----
  # One per element of the coefficient vector, not one per selected position:
  # the sampler masks the regressors with diag(lambda) as a whole, so a position
  # the sweep never visits keeps whatever it starts with for the entire run.
  # Starting every position at one therefore leaves the unselected coefficients
  # in the model, which is what they are. That covers the k * r loadings at the
  # front of 'a' in particular, which a VEC never selects over --
  # inclusion_prior() drops them from 'include' and the sampler refuses them
  # there -- and which a zero here would switch off permanently.
  use_varsel <- object[["model"]][["varsel"]] %in% c("ssvs", "bvs")

  if (use_varsel & !is.null(object[["data"]][["train"]][["z"]])) {
    object[["initial"]][["a_lambda"]] <- matrix(1, ncol(object[["data"]][["train"]][["z"]]))
  }
  if (use_varsel & !is.null(object[["priors"]][["psi"]][["inprior"]])) {
    object[["initial"]][["psi_lambda"]] <- matrix(1, nrow(object[["priors"]][["psi"]][["inprior"]]))
  }

  # Variances of state equations ----
  object <- .add_initial_values_state_errors(object, method = method)
  
  # Initial values for errors
  object <- .add_initial_values_measurement_errors(object = object,
                                                   method = method,
                                                   u = u)
  
  return(object)
}



# Square root of a matrix
# Used in add_initial_values methods to estimate VEC models
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



# Johansen's maximum likelihood estimate of the error correction term
#
# Used for the initial values of a VEC model and for a cointegration space prior
# centred on it (add_priors with coint$p_tau_i = "ml"). Works on
# object$data$train as it stands, so on the scaled series if
# scale_error_correction() has been applied.
#
# Returns the M x r cointegration matrix 'beta', normalised to beta' S11 beta =
# I, the k x r loadings 'alpha' belonging to it, the k x k residual covariance
# 'omega', and the M x T residuals 'r1' of the regression of w on the
# short-run regressors, whose cross product is the information the estimate of
# beta carries.
.coint_ml <- function(object) {

  y <- t(object[["data"]][["train"]][["y"]])
  w <- t(object[["data"]][["train"]][["w"]])
  x <- object[["data"]][["train"]][["x"]]
  r <- object[["model"]][["rank"]]
  tt <- ncol(y)

  if (!is.null(x) && NCOL(x) > 0) {
    x <- t(x)
    M <- diag(tt) - crossprod(x, solve(tcrossprod(x))) %*% x
    R0 <- y %*% M # Residuals of regression of y on x
    R1 <- w %*% M # Residuals of regression of w on x
  } else {
    R0 <- y
    R1 <- w
  }
  S00 <- tcrossprod(R0) / tt
  S00_inv <- solve(S00)
  S01 <- tcrossprod(R0, R1) / tt
  S10 <- tcrossprod(R1, R0) / tt
  S11 <- tcrossprod(R1) / tt
  S11_sqrt_inv <- solve(.mroot(S11))
  lambda <- eigen(S11_sqrt_inv %*% S10 %*% S00_inv %*% S01 %*% t(S11_sqrt_inv))

  beta <- t(crossprod(matrix(lambda$vectors[, 1:r], nrow(w)), S11_sqrt_inv))
  # beta' S11 beta is the identity up to rounding; not relied upon here.
  bsb <- crossprod(beta, S11 %*% beta)
  alpha <- S01 %*% beta %*% solve(bsb)
  omega <- S00 - alpha %*% bsb %*% t(alpha)

  return(list(beta = beta, alpha = alpha, omega = omega, r1 = R1))
}
