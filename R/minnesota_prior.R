#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a VAR model.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' or \code{\link{gen_vec}}.
#' @param kappa0 a numeric specifying the prior standard deviation of coefficients that correspond to
#' own lags of endogenous variables.
#' @param kappa1 a numeric specifying the size of the prior standard deviations of endogenous
#' variables, which do not correspond to own lags, relative to argument \code{kappa0}.
#' @param kappa2 a numeric specifying the size of the prior standard deviations of exogenous
#' variables relative to argument \code{kappa0}.
#' @param kappa3 a numeric specifying the size of the prior standard deviations of deterministic
#' terms relative to argument \code{kappa0}.
#' @param max_var a positive numeric specifying the maximum prior variance that is allowed for
#' coefficients of non-deterministic variables. If \code{NULL} (default), the prior variances are not limited.
#' @param coint_var a logical specifying whether the model is a cointegrated VAR model,
#' for which the prior means of first own lags should be set to one.
#' 
#' @details The function calculates the Minnesota prior of a VAR model. For the endogenous variable
#' \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#' \deqn{ \left( \frac{\kappa_{0}}{l} \right)^2 \ \text{for own lags of endogenous variables,}}
# \deqn{ \left( \frac{\kappa_0 \kappa_1}{l} \frac{\sigma_{i}}{\sigma_{j}} \right)^2 \text{ for endogenous variables other than own lags,}}
# \deqn{ \left( \frac{\kappa_0 \kappa_2}{l} \frac{\sigma_{i}}{\sigma_{j}} \right)^2  \text{ for exogenous variables,}}
# \deqn{ (\kappa_0 \kappa_3)^2 \text{ for deterministic terms,}}
#' where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#' OLS estimate of the model. For exogenous variables \eqn{\sigma_{i}} corresponds to the standard
#' deviation of the original series.
#' 
#' For VEC models the function only provides priors for the non-cointegration part of the model. The
#' residual standard errors \eqn{\sigma_i} are based on an unrestricted OLS regression of the
#' endogenous variables on the error correction term and the non-cointegration regressors.
#' 
#' @return A list containing a matrix of prior means and the precision matrix.
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
minnesota_prior <- function(object, kappa0 = 2, kappa1 = .5, kappa2 = .5, kappa3 = 5, max_var = NULL, coint_var = FALSE) {
  
  if (any(c(kappa0, kappa1, kappa2, kappa3) <= 0)) {
    "Kappa arguments must be positive."
  }
  
  if (!is.null(max_var)) {
    if (max_var <= 0) {
      "Argument 'max_var' must be positive."
    } 
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
  
  V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
  
  # Endogenous variables
  s_endo <- y %*% (diag(1, tt) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / tt
  s_endp <- sqrt(diag(s_endo)) # Residual standard deviations (OLS)
  if (p > 0) {
    for (r in 1:p) {
      for (l in 1:k) {
        for (j in 1:k) {
          if (l == j) {
            V[l, n_ect + (r - 1) * k + j] <- (kappa0 / r)^2
          } else {
            V[l, n_ect + (r - 1) * k + j] <- (kappa0 * kappa1 / r * s_endo[l] / s_endo[j])^2
          }
        } 
      }
    } 
  }
  
  # Exogenous variables
  m <- 0
  s <- 0
  if (!is.null(object$model$exogen)) {
    m <- length(object$model$exogen$variables)
    s <- object$model$exogen$lags
    if (type == "VEC") {
      s <- s - 1
    }
    s_exo <- sqrt(apply(x[p * k + 1:m,], 1, stats::var))
    for (r in 1:s) {
      for (l in 1:k) {
        for (j in 1:m) {
          # Note that in the loop r starts at 1, so that this is equivalent to l + 1
          V[l, n_ect + p * k + (r - 1) * m + j] <- (kappa0 * kappa2 / r * s_endo[l] / s_exo[j])^2
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
    V[, -(1:(n_ect + k * p + m * s))] <- (kappa0 * kappa3)^2
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
                 "v_i" = v_i)
  
  return(result)
}