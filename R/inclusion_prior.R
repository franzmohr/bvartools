#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' or \code{\link{gen_vec}}.
#' @param prob a numeric specifying the prior inclusion probability of all model parameters.
#' @param minnesota_like logical. If \code{TRUE}, the prior inclusion probabilities of the
#' parameters are calculated in a similar way as the Minnesota prior. See 'Details'.
#' @param kappa a numeric vector of four elements containing the prior inclusion probabilities
#' of coefficients that correspond to own lags of endogenous variables, to endogenous variables,
#' which do not correspond to own lags, to exogenous variables and deterministic terms, respectively.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param exclude_deterministics logical. If \code{TRUE} (default), the vector of variables, which should
#' be included in the BVS algorithm does not include deterministic terms.
#' 
#' @details Prior inclusion probabilities \eqn{\underline{\pi}_1} are calculated as
#' \tabular{c l}{
#' \eqn{\frac{\kappa_1}{r}} & \text{for own lags of endogenous variables,} \\
#' \eqn{\frac{\kappa_2}{r}} & \text{for other endogenous variables,} \\
#' \eqn{\frac{\kappa_3}{1 + r}} & \text{for exogenous variables,}\\
#' \eqn{\kappa_{4}} & \text{for deterministic variables,} 
#' \end{align}}
#' for lag \eqn{r} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' thrid and forth element in \code{kappa}, respectively.
#' 
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selction algorithm.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \url{https://doi.org/10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \url{https://doi.org/10.1002/jae.1271}
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
#' # Obtain inclusion prior
#' pi_prior <- inclusion_prior(object)
#' 
#' @export
inclusion_prior <- function(object, prob = .5, minnesota_like = FALSE,
                      kappa = c(0.8, 0.5, 0.5, .8), exclude_deterministics = TRUE) {
  
  if (prob > 1 | prob < 0) {
    stop("Argument 'prob' must be between 0 and 1.")
  }
  if (any(kappa > 1) | any(kappa < 0)) {
    stop("Argument 'kappa' may only contain values between 0 and 1.")
  }
  
  y <- object$Y
  tt <- NCOL(y)
  
  k <- length(object$model$endogenous$variables)
  m <- 0
  n <- 0
  p <- object$model$endogenous$lags
  s <- 0
  
  if (!is.null(object$model$exogen)) {
    m <- length(object$model$exogen$variables)
    s <- object$model$exogen$lags
  }

  if (object$model$type == "VAR") {
    x <- object$Z 
    if (!is.null(object$model$deterministic)) {
      n <- length(object$model$deterministic)
    }
  }
  
  if (object$model$type == "VEC") {
    if (!is.null(object$X)) {
      x <- rbind(object$W, object$X)
    } else {
      x <- object$W
    }
    if (!is.null(object$model$deterministic$unrestricted)) {
      n <- n + length(object$model$deterministic$unrestricted)
    }
    p <- p - 1
    s <- s - 1
  }
  
  result <- matrix(NA, k, NROW(x))
  include <- matrix(1:(k * NROW(x)))
  
  if (minnesota_like) {
    if (p > 0) {
      for (i in 1:p) {
        result[, (i - 1) * k + 1:k] <- kappa[2] / i
        if (k > 1) {
          diag(result[, (i - 1) * k + 1:k]) <- kappa[1] / i 
        } else {
          result[1, (i - 1) * k + 1:k] <- kappa[1] / i
        }
      }
    }
    
    if (m > 0) {
      result[, p * k + 1:m] <- kappa[3]
      if (s > 0) {
        for (i in 1:s) {
          result[, p * k + m + (i - 1) * k + 1:k] <- kappa[3] / (1 + i)
        }
      }
    }
    
    if (n > 0) {
      result[, p * k + (s + 1) * m + 1:n] <- kappa[4]
    }
  } else {
    result[,] <- prob
  }
  
  result <- matrix(result)
  
  if (n > 0 & exclude_deterministics) {
    pos_det <- k * (k * p + (s + 1) * m) + 1:(k * n)
    include <- matrix(include[-pos_det, ])
    if (NROW(include) == 0) {
      warning("Deterministic terms are excluded from BVS and no further parameters are estimated.")
    }
  }
  
  result <- list("prior" = result,
                 "include" = include)
  
  return(result)
}