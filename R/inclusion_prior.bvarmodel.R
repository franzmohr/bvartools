#' Prior Inclusion Probabilities
#' 
#' Prior inclusion probabilities as required for stochastic search variable selection (SSVS) à la
#' George et al. (2008) and Bayesian variable selection (BVS) à la Korobilis (2013).
#' 
#' @param object an object of class 'bvarmodel', usually, a result of a
#' call to \code{\link{create_bvarmodel}}.
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
#' \eqn{\frac{\kappa_3}{1 + r}} \tab for unmodelled exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{r} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' third and forth element in \code{kappa}, respectively.
#' 
#' @return A list containing a matrix of prior inclusion probabilities and an integer vector
#' specifying the positions of variables, which should be included in the variable selection algorithm.
#' 
#' @examples
#' 
#' # Prepare data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model input
#' object <- create_bvarmodel(e1)
#' 
#' # Obtain inclusion prior
#' incl <- inclusion_prior(object)
#' 
#' @export
inclusion_prior.bvarmodel <- function(object, prob = .5, exclude_deterministics = TRUE,
                                      minnesota_like = FALSE, kappa = c(0.8, 0.5, 0.5, 0.8)) {
  
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
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    
    z <- object[["data"]][["train"]][["z"]]
    k <- object[["model"]][["k"]]
    tt <- nrow(object[["data"]][["train"]][["y"]])
    p <- object[["model"]][["p"]]
    n_a <- k * p
    m <- object[["model"]][["m"]]
    s <- object[["model"]][["s"]]
    n_b <- m * (s + 1)
    n_c <- object[["model"]][["n"]]
    
    inprior <- rep(prob, ncol(z))
    include <- 1:ncol(z)
    
    if (minnesota_like) {
      
      incl_matrix <- matrix(NA, k, n_a + n_b + n_c)
      
      if (p > 0) {
        for (i in 1:p) {
          incl_matrix[, (i - 1) * k + 1:k] <- kappa[2] / i
          if (k > 1) {
            diag(incl_matrix[, (i - 1) * k + 1:k]) <- kappa[1] / i 
          } else {
            incl_matrix[, (i - 1) * k + 1] <- kappa[1] / i
          }
        }
      }
      
      if (m > 0) {
        incl_matrix[, n_a + 1:m] <- kappa[3]
        if (s > 0) {
          for (i in 1:s) {
            incl_matrix[, n_a + m + (i - 1) * m + 1:m] <- kappa[3] / (1 + i)
          }
        }
      }
      
      if (n_c > 0) {
        incl_matrix[, n_a + n_b + 1:n_c] <- kappa[4]
      }
      
      inprior[1:(k * (n_a + n_b + n_c))] <- c(incl_matrix)
    }
    
    # Exclude deterministics from variables selection algorithm
    if (n_c > 0 & exclude_deterministics) {
      pos_det <- k * (n_a + n_b) + 1:(k * n_c)
      include <- include[-pos_det]
    }
    
    if (length(include) > 0) {
      result <- list("prior" = matrix(inprior),
                     "include" = matrix(include))
    }
  }
  
  return(result)
}