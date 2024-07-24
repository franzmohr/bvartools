#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class \code{"bvecmodel"}, usually, a result of a
#' call to \code{\link{create_vec_model}}.
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
#' specifying the positions of variables, which should be included in the variable selection algorithm.
#' 
#' @examples
#' 
#' # Prepare data
#' data("e6")
#' 
#' # Generate model input
#' object <- create_vec_model(e6, r = 1)
#' 
#' # Obtain inclusion prior
#' pi_prior <- inclusion_prior(object)
#' 
#' @export
inclusion_prior.bvecmodel <- function(object, prob = .5, exclude_deterministics = TRUE,
                                      minnesota_like = FALSE, kappa = c(0.8, 0.5, 0.5, .8)) {
  
  # Input checks
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
  
  result <- NULL
  if (!is.null(object[["data"]][["z"]])) {
    
    z <- object[["data"]][["z"]]
    k <- object$model$k
    r <- object$model$rank
    n_c_unres <- object$model$n_unrestricted
    n_alpha <- k * r
    
    inprior <- rep(prob, ncol(z))
    exclude <- NULL
    include <- 1:ncol(z)
    
    # alpha coefficients ----
    if (object$model$rank > 0) {
      inprior[1:n_alpha] <- NA
      exclude <- append(exclude, 1:n_alpha)
    }
    
    # non-alpha coefficients ----
    if (minnesota_like & !is.null(object$data$x)) {
      
      p <- object[["model"]][["p"]]
      n_gamma <- k * (p - 1)
      m <- object[["model"]][["m"]]
      s <- object[["model"]][["s"]]
      n_upsilon <- m * s
      
      incl_matrix <- matrix(NA, k, n_gamma + n_upsilon + n_c_unres)
      
      if (p > 1) {
        for (i in 1:(p - 1)) {
          incl_matrix[, (i - 1) * k + 1:k] <- kappa[2] / i
          if (k > 1) {
            diag(incl_matrix[, (i - 1) * k + 1:k]) <- kappa[1] / i 
          } else {
            incl_matrix[, (i - 1) * k + 1] <- kappa[1] / i
          }
        }
      }
      
      if (m > 0) {
        incl_matrix[, n_gamma + 1:m] <- kappa[3]
        if (s > 0) {
          for (i in 1:s) {
            incl_matrix[, n_gamma + m + (i - 1) * m + 1:m] <- kappa[3] / (1 + i)
          }
        }
      }
      
      if (n_c_unres > 0) {
        incl_matrix[, n_gamma + n_upsilon + 1:n_c_unres] <- kappa[4]
      }
      
      inprior[n_alpha + 1:(k * (n_gamma + n_upsilon + n_c_unres))] <- c(incl_matrix)
    } 
    
    # Exclude deterministics from variables selection algorithm
    if (n_c_unres > 0 & exclude_deterministics) {
      exclude <- append(exclude, n_alpha + k * (n_gamma + n_upsilon) + 1:(k * n_c_unres))
      include <- include[-exclude]
    }
    
    if (length(include) > 0) {
      result <- list("prior" = matrix(inprior),
                     "include" = matrix(include))
    }
  }
  
  return(result)
}
