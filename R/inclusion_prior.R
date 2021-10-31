#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' or \code{\link{gen_vec}}.
#' @param prob a numeric specifying the prior inclusion probability of all model parameters.
#' @param exclude_deterministics logical. If \code{TRUE} (default), the vector of the positions of
#' included variables does not include the positions of deterministic terms.
#' @param minnesota_like logical. If \code{TRUE}, the prior inclusion probabilities of the
#' parameters are calculated in a similar way as the Minnesota prior. See 'Details'.
#' @param kappa a numeric vector of four elements containing the prior inclusion probabilities
#' of coefficients that correspond to own lags of endogenous variables, to endogenous variables,
#' which do not correspond to own lags, to exogenous variables and deterministic terms, respectively.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' 
#' @details If \code{minnesota_like = TRUE}, prior inclusion probabilities \eqn{\underline{\pi}_1}
#' are calculated as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{r}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{r}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + r}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{r} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' third and forth element in \code{kappa}, respectively.
#' 
#' For vector error correction models the function generates prior inclusion probabilities for differenced
#' variables and unrestricted deterministc terms as described above. For variables in the error
#' correction term prior inclusion probabilites are calculated as
#' \tabular{cl}{
#' \eqn{\kappa_1} \tab fow own levels of endogenous variables, \cr
#' \eqn{\kappa_2} \tab for levels of other endogenous variables, \cr
#' \eqn{\kappa_3} \tab for levels of exogenous variables, \cr
#' \eqn{\kappa_4} \tab for deterministic variables.
#' }
#' 
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selction algorithm.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Prepare data
#' data("e1")
#' 
#' # Generate model input
#' object <- gen_var(e1)
#' 
#' # Obtain inclusion prior
#' pi_prior <- inclusion_prior(object)
#' 
#' @export
inclusion_prior <- function(object, prob = .5, exclude_deterministics = TRUE,
                            minnesota_like = FALSE, kappa = c(0.8, 0.5, 0.5, .8)) {
  
  if (!minnesota_like) {
    if (prob > 1 | prob < 0) {
      stop("Argument 'prob' must be between 0 and 1.")
    } 
  }
  if (minnesota_like) {
    if (any(kappa > 1) | any(kappa < 0)) {
      stop("Argument 'kappa' may only contain values between 0 and 1.")
    } 
  }
  
  y <- t(object$data$Y)
  tt <- NCOL(y)
  
  k <- length(object$model$endogen$variables)
  m <- 0
  n <- 0
  n_res <- 0
  p <- object$model$endogen$lags
  s <- 0
  n_pi <- 0
  
  if (!is.null(object$model$exogen)) {
    m <- length(object$model$exogen$variables)
    s <- object$model$exogen$lags
  }

  if (object$model$type == "VAR") {
    x <- t(object$data$Z)
    tot_par <- ncol(object[["data"]][["SUR"]])
    if (!is.null(object$model$deterministic)) {
      n <- length(object$model$deterministic)
    }
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
    tot_par <- ncol(object[["data"]][["SUR"]])
    
    if (object$model$rank != 0) {
      if (!is.null(object$model$deterministic$restricted)) {
        n_res <- length(object$model$deterministic$restricted)
      } 
    }
    if (!is.null(object$model$deterministic$unrestricted)) {
      n <- n + length(object$model$deterministic$unrestricted)
    }
    p <- p - 1
    s <- s - 1
    if (object$model$rank != 0) {
      n_pi <- k + m + n_res
    }
  }
  
  result <- matrix(NA, k, NROW(x))
  include <- matrix(1:(k * NROW(x)))
  
  if (minnesota_like) {
    if (object$model$type == "VEC") {
      if (object$model$rank == 0) {
        result[, 1:k] <- kappa[2]
        if (k > 1) {
          diag(result[, 1:k]) <- kappa[1] 
        } else {
          result[, 1] <- kappa[1]
        }
        if (m > 0) {
          result[, k + 1:m] <- kappa[3]
        }
        if (n_res > 0) {
          result[, k + m + 1:n_res] <- kappa[4]
        } 
      }
    }
    if (p > 0) {
      for (i in 1:p) {
        result[, n_pi + (i - 1) * k + 1:k] <- kappa[2] / i
        if (k > 1) {
          diag(result[, n_pi + (i - 1) * k + 1:k]) <- kappa[1] / i 
        } else {
          result[, n_pi + (i - 1) * k + 1] <- kappa[1] / i
        }
      }
    }
    
    if (m > 0) {
      result[, n_pi + p * k + 1:m] <- kappa[3]
      if (s > 0) {
        for (i in 1:s) {
          result[, n_pi + p * k + m + (i - 1) * m + 1:m] <- kappa[3] / (1 + i)
        }
      }
    }
    
    if (n > 0) {
      result[, n_pi + p * k + (s + 1) * m + 1:n] <- kappa[4]
    }
  } else {
    result[,] <- prob
  }
  
  result <- matrix(result)
  if (object[["model"]][["structural"]]) {
    result <- rbind(result, matrix(prob, k * (k - 1) / 2))
    include <- rbind(include, matrix(k * NROW(x) + 1:(k * (k - 1) / 2)))
  }
  
  # Exclude deterministics from variables selection algorithm
  if ((n > 0 | n_res > 0) & exclude_deterministics) {
    if (object$model$type == "VAR") {
      pos_det <- k * (k * p + (s + 1) * m) + 1:(k * n) 
    }
    if (object$model$type == "VEC") {
      pos_det <- NULL
      if (n_res > 0 & n_pi > 0) {
        pos_det <- c(pos_det, k * (k + m) + 1:(k * n_res))
      }
      if (n > 0) {
        pos_det <- c(pos_det, k * (n_pi + k * p + m * (s + 1)) + 1:(k * n))
      }
    }
    include <- matrix(include[-pos_det, ])
    if (NROW(include) == 0) {
      warning("Deterministic terms are excluded from the variable selection algorithm and no further parameters are estimated.")
    }
  }
  
  result <- list("prior" = result,
                 "include" = include)
  
  return(result)
}