#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a VAR model.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' or \code{\link{gen_vec}}.
#' @param kappa0 a numeric specifying the prior variance of coefficients that correspond to
#' own lags of endogenous variables.
#' @param kappa1 a numeric specifying the size of the prior variance of endogenous
#' variables, which do not correspond to own lags, relative to argument \code{kappa0}.
#' @param kappa2 a numeric specifying the size of the prior variance of non-deterministic exogenous
#' variables relative to argument \code{kappa0}. Default is \code{NULL}, which indicates that the formula
#' for the calculation of the prior variance of deterministic terms is used for all exogenous variables.
#' @param kappa3 a numeric specifying the size of the prior variance of deterministic
#' terms relative to argument \code{kappa0}.
#' @param max_var a positive numeric specifying the maximum prior variance that is allowed for
#' coefficients of non-deterministic variables. If \code{NULL} (default), the prior variances are not limited.
#' @param coint_var a logical specifying whether the model is a cointegrated VAR model,
#' for which the prior means of first own lags should be set to one.
#' @param sigma either \code{"AR"} (default) or \code{"VAR"} indicating that the variances of the endogenous
#' variables \eqn{\sigma^2} are calculated based on a univariate AR regression or a least squares estimate of
#' the VAR form, respectively. In both cases all deterministic variables are used in the regressions,
#' if they appear in the model.
#' 
#' @details The function calculates the Minnesota prior of a VAR model. For the endogenous variable
#' \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#' \deqn{ \frac{\kappa_{0}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#' \deqn{ \frac{\kappa_{0} \kappa_{1}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#' \deqn{ \frac{\kappa_{0} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#' \deqn{ \kappa_{0} \kappa_{3} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#' where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#' LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#' 
#' For VEC models the function only provides priors for the non-cointegration part of the model. The
#' residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#' endogenous variables on the error correction term and the non-cointegration regressors.
#' 
#' @return A list containing a matrix of prior means and the precision matrix of the cofficients and the
#' inverse variance-covariance matrix of the error term.
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
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
#' # Obtain Minnesota prior
#' prior <- minnesota_prior(object)
#' 
#' @export
minnesota_prior <- function(object, kappa0 = 2, kappa1 = .5, kappa2 = NULL, kappa3 = 5,
                            max_var = NULL, coint_var = FALSE, sigma = "AR") {
  
  if (any(c(kappa0, kappa1, kappa2, kappa3) <= 0)) {
    stop("Kappa arguments must be positive.")
  }
  
  if (!is.null(max_var)) {
    if (max_var <= 0) {
      stop("Argument 'max_var' must be positive.")
    } 
  }
  
  if (!sigma %in% c("AR", "VAR")) {
    stop("Argument 'sigma' must be either 'AR' or 'VAR'.")
  }
  
  y <- object$Y
  type <- object$model$type
  if (type == "VAR") {
    x <- object$Z 
  }
  
  n_ect <- 0
  if (type == "VEC") {
    if (!is.null(object$X)) {
      x <- rbind(object$W, object$X)
    } else {
      x <- object$W 
    }
    n_ect <- NROW(object$W)
  }
  k <- NROW(y)
  tt <- NCOL(y)
  tot_par <- k * NROW(x)
  p <- object$model$endogenous$lags
  if (type == "VEC") {
    p <- p - 1
  }
  
  m <- 0
  s <- 0
  if (!is.null(object$model$exogen)) {
    m <- length(object$model$exogen$variables)
    s <- object$model$exogen$lags
  }
  
  V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
  
  # Obtain OLS sigma
  ols_sigma <- y %*% (diag(1, tt) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / (tt - nrow(x))
  
  # Determine positions of deterministic terms for calculation of sigma
  pos_det <- NULL
  if (!is.null(object$model$deterministic)) {
    if (type == "VAR") {
      pos_det <- k * p + s * m + 1:length(object$model$deterministic) 
    }
    if (type == "VEC") {
      if (!is.null(object$model$deterministic$restricted)) {
        pos_det <- c(pos_det, k + m + 1:length(object$model$deterministic$restricted))
      }
      if (!is.null(object$model$deterministic$unrestricted)) {
        pos_det <- c(pos_det, n_ect + k * p + m * s + 1:length(object$model$deterministic$unrestricted))
      }
    }
  }
  
  # Obtain sigmas for V_i
  if (sigma == "AR") { # Univariate AR
    s_endo <- diag(0, k)
    for (i in 1:k) {
      if (type == "VAR") {
        pos <- c(i + k * ((1:p) - 1), pos_det)
      }
      if (type == "VEC") {
        if (p > 0) {
          pos <- c(i, n_ect + i + k * ((1:p) - 1), pos_det)
        } else {
          pos <- c(i, pos_det)
        }
      }
      y_temp <- matrix(y[i, ], 1)
      x_temp <- matrix(x[pos, ], length(pos))
      s_endo[i, i] <- y_temp %*% (diag(1, tt) - t(x_temp) %*% solve(tcrossprod(x_temp)) %*% x_temp) %*% t(y_temp) / (tt - length(pos))
    }
  }
  if (sigma == "VAR") { # VAR model
    s_endo <- ols_sigma
  }
  s_endo <- sqrt(diag(s_endo)) # Residual standard deviations (OLS)
  
  # Endogenous variables
  if (p > 0) {
    for (r in 1:p) {
      for (l in 1:k) {
        for (j in 1:k) {
          if (l == j) {
            V[l, n_ect + (r - 1) * k + j] <- kappa0 / r^2
          } else {
            V[l, n_ect + (r - 1) * k + j] <- kappa0 * kappa1 / r^2 * s_endo[l]^2 / s_endo[j]^2
          }
        } 
      }
    } 
  }
  
  # Exogenous variables
  if (!is.null(object$model$exogen)) {
    if (type == "VEC") {
      s <- s - 1
    }
    s_exo <- sqrt(apply(x[p * k + 1:m,], 1, stats::var))
    for (r in 1:s) {
      for (l in 1:k) {
        for (j in 1:m) {
          # Note that in the loop r starts at 1, so that this is equivalent to l + 1
          if (!is.null(kappa2)) {
            V[l, n_ect + p * k + (r - 1) * m + j] <- kappa0 * kappa2 / r^2 * s_endo[l]^2 / s_exo[j]^2 
          } else {
            V[l, n_ect + p * k + (r - 1) * m + j] <- kappa0 * kappa3 * s_endo^2
          }
        }
      }
    } 
  }
  
  # Restrict prior variances
  if (!is.null(max_var)) {
    if (any(V > max_var)) {
      V[which(V > max_var)] <- max_var
    } 
  }
  
  # Deterministic variables
  if (!is.null(object$model$deterministic)){
    for (i in 1:k) {
      V[, -(1:(n_ect + k * p + m * s))] <- kappa0 * kappa3 * s_endo^2 
    }
  }
  
  # Drop cointegration priors
  if (type == "VEC") {
    V <- V[, -(1:n_ect)]
    tot_par <- k * NROW(object$X)
  }
  
  # Prior means
  mu <- matrix(rep(0, tot_par), k)
  if (coint_var & object$model$type == "VAR") {
    mu[1:k, 1:k] <- diag(1, k)
  }
  mu <- matrix(mu)
  
  # Prior precision
  v_i <- diag(c(1 / V))
  
  result <- list("mu" = mu,
                 "v_i" = v_i,
                 "sigma_i" = solve(ols_sigma))
  
  return(result)
}