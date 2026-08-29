#' Minnesota Prior
#' 
#' Calculates the Minnesota prior for a VEC model.
#' 
#' @param object an object of class 'bvecmodel', usually, a result of a call to
#' \code{\link{create_bvecmodel}}.
#' @param kappa1 a numeric specifying the prior variance of coefficients that correspond to
#' own lags of endogenous variables. See 'Details'.
#' @param kappa2 a numeric specifying the size of the prior variance of endogenous
#' variables, which do not correspond to own lags. See 'Details'.
#' @param kappa3 a numeric specifying the size of the prior variance of non-deterministic exogenous
#' variables. Default is \code{NULL}, which indicates that the formula
#' for the calculation of the prior variance of deterministic terms is used
#' for all exogenous variables. See 'Details'.
#' @param kappa4 a numeric specifying the size of the prior variance of deterministic
#' terms. See 'Details'.
#' @param max_var a positive numeric specifying the maximum prior variance that is allowed for
#' coefficients of non-deterministic variables. If \code{NULL} (default), the prior variances are not limited.
#' @param sigma either \code{"AR"} (default) or \code{"VAR"} indicating that the variances of the endogenous
#' variables \eqn{\sigma^2} are calculated based on a univariate AR regression or a least squares estimate of
#' the VAR form, respectively. In both cases all deterministic variables are used in the regressions,
#' if they appear in the model.
#' 
#' @details The function calculates the Minnesota prior in the same manner as for
#' a traditional VAR model. For the endogenous variable
#' \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#' \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#' \deqn{ \frac{\kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#' \deqn{ \frac{\kappa_{3}}{(l + 1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#' \deqn{ \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#' where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#' LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#' 
#' The function only provides priors for the non-cointegration part of the model. The
#' residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#' endogenous variables on the error correction term and the non-cointegration regressors.
#' 
#' @return A list containing a matrix of prior means and the precision matrix of the coefficients and the
#' inverse variance-covariance matrix of the error term, which was obtained by an LS estimation.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2020). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples
#' 
#' # Load data
#' data("e6")
#' 
#' # Generate model input
#' object <- create_bvecmodel(e6, r = 1)
#' 
#' # Obtain Minnesota prior
#' prior <- minnesota_prior(object)
#' 
#' @export
minnesota_prior.bvecmodel <- function(object, kappa1 = 2, kappa2 = 1, kappa3 = NULL, kappa4 = 10,
                                      max_var = NULL, sigma = "AR") {
  
  if (any(c(kappa1, kappa2, kappa3, kappa4) <= 0)) {
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
  
  y <- t(object$data$y)
  type <- object$model$type
  k <- NROW(y)

  mu <- NULL  
  V <- NULL
  result <- NULL
  
  if (!is.null(object$data$z)) {
    
    if (!is.null(object$data$x)) {
      
      
      x <- t(cbind(object$data$w, object$data$x))
      n_ect <- NCOL(object$data$w)
      tt <- NCOL(y)
      tot_par <- k * NROW(x)
      p <- object$model$p
      p <- p - 1
      
      m <- 0
      s <- 0
      if (object$model$m > 0) {
        m <- object$model$m
        s <- object$model$s
      }
      
      V <- matrix(rep(NA, tot_par), k) # Set up matrix for variances
      
      # Obtain OLS sigma
      ols_sigma <- y %*% (diag(1, tt) - t(x) %*% solve(tcrossprod(x)) %*% x) %*% t(y) / (tt - nrow(x))
      
      # Determine positions of deterministic terms for calculation of sigma
      pos_det <- NULL
      if (object$model$n_restricted > 0 | object$model$n_unrestricted > 0) {
        if (object$model$n_restricted > 0) {
          pos_det <- c(pos_det, k + m + 1:length(object$model$n_restricted))
        }
        if (object$model$n_unrestricted > 0) {
          pos_det <- c(pos_det, n_ect + k * p + m * s + 1:length(object$model$n_unrestricted))
        }
      }
      
      # Obtain sigmas for V_i
      if (sigma == "AR") { # Univariate AR
        s_endo <- diag(0, k)
        if (p > 0 | !is.null(pos_det)) {
          for (i in 1:k) {
            
            if (p > 0) {
              pos <- c(i, n_ect + i + k * ((1:p) - 1), pos_det)
            } else {
              pos <- c(i, pos_det)
            }
            
            y_temp <- matrix(y[i, ], 1)
            x_temp <- matrix(x[pos,], length(pos))
            s_endo[i, i] <- y_temp %*% (diag(1, tt) - t(x_temp) %*% solve(tcrossprod(x_temp)) %*% x_temp) %*% t(y_temp) / (tt - length(pos))
          } 
        } else {
          diag(s_endo) <- apply(matrix(y, k), 1, stats::var)
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
                V[l, n_ect + (r - 1) * k + j] <- kappa1 / r^2
              } else {
                V[l, n_ect + (r - 1) * k + j] <- kappa2 / r^2 * s_endo[l]^2 / s_endo[j]^2
              }
            } 
          }
        } 
      }
      
      # Exogenous variables
      if (m > 0) {
        s <- s - 1
        s_exo <- sqrt(apply(matrix(x[p * k + 1:m,], m), 1, stats::var))
        for (r in 1:(s + 1)) {
          for (l in 1:k) {
            for (j in 1:m) {
              # Note that in the loop r starts at 1, so that this is equivalent to l + 1
              if (is.null(kappa3)) {
                V[l, n_ect + p * k + (r - 1) * m + j] <- kappa4 * s_endo[l]^2
              } else {
                V[l, n_ect + p * k + (r - 1) * m + j] <- kappa3 / r^2 * s_endo[l]^2 / s_exo[j]^2 
              }
            }
          }
        } 
      }
      
      # Restrict prior variances
      if (!is.null(max_var)) {
        if (any(stats::na.omit(c(V)) > max_var)) {
          V[which(V > max_var)] <- max_var
        } 
      }
      
      # Deterministic variables
      if (object$model$n_unrestricted > 0){
        for (i in 1:k) {
          V[, -(1:(n_ect + k * p + m * s))] <- kappa4 * s_endo^2 
        }
      }
      
      # Drop cointegration priors
      V <- V[, -(1:n_ect)]
      tot_par <- k * NCOL(object$data$x)
      
      
      # Prior means
      mu <- matrix(rep(0, tot_par), k)
      mu <- matrix(mu)
      
      V <- matrix(V)
      
    }else {
      s_endo <- sqrt(diag(stats::var(t(y))))
    }
    
    
    # Structural parameters
    if (object$model$structural & k > 1) {
      mu <- rbind(mu, matrix(0, k * (k - 1) / 2))
      
      V_struct <- matrix(NA, k, k)
      for (j in 1:(k - 1)) {
        V_struct[(j + 1):k, j] <- kappa2 * s_endo[(j + 1):k]^2 / s_endo[j]^2  
      }
      V_struct <- matrix(V_struct[lower.tri(V_struct)])
      V <- rbind(V, V_struct)
    }
    
    if (object[["model"]][["rank"]] > 0) {
      n_alpha <- k * object[["model"]][["rank"]]
      mu <- rbind(mu, matrix(0, n_alpha))
      V <- rbind(V, matrix(kappa0, n_alpha))
    }
    
    # Prior precision
    v_i <- diag(c(1 / V))
    
    result <- list("mu" = mu,
                   "v_i" = v_i)
    
    if (!is.null(object[["data"]][["train"]][["x"]])) {
      if (sigma == "AR") {
        result[["sigma_inv"]] <- matrix(0, k, k)
        diag(result[["sigma_inv"]]) = 1 / s_endo^2
      }
      if (sigma == "VAR") {
        result[["sigma_inv"]] = solve(ols_sigma)  
      }
    }
    
    if (length(result$mu) == 0) {
      result <- NULL
    }
    
  }
  
  return(result)
}