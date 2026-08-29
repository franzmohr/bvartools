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
#' @param kappa1 a numeric specifying the prior inclusion probability of
#' coefficients that correspond to own lags of endogenous variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa2 a numeric specifying the size of the prior inclusion probabilities
#' of endogenous variables, which do not correspond to own lags.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa3 a numeric specifying the size of the prior inclusion probabilities
#' of non-deterministic exogenous variables. Default is \code{NULL}, which indicates that the formula
#' for the calculation of the prior inclusion probabilities of deterministic terms
#' is used for all exogenous variables.
#' Only used if \code{minnesota_like = TRUE}. See 'Details'.
#' @param kappa4 a numeric specifying the size of the prior inclusion probabilities
#' of deterministic terms. Only used if \code{minnesota_like = TRUE}. See 'Details'.
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
inclusion_prior.bvarmodel <- function(object,
                                      prob = .5,
                                      exclude_deterministics = TRUE,
                                      minnesota_like = FALSE,
                                      kappa1 = 0.8,
                                      kappa2 = 0.5,
                                      kappa3 = 0.5,
                                      kappa4 = 0.8) {
  
  if (!minnesota_like) {
    if (prob > 1 | prob < 0) {
      stop("Argument 'prob' must be between 0 and 1.")
    } 
  }
  if (minnesota_like) {
    
    if (kappa1 < 0) {
      stop("Argument 'kappa1' must not be negative.")
    }
    if (kappa1 > 1) {
      stop("Argument 'kappa1' must not be larger than 1.")
    }
    
    if (kappa2 <= 0) {
      stop("Argument 'kappa2' must not be negative.")
    }
    if (kappa2 > 1) {
      stop("Argument 'kappa2' must not be larger than 1.")
    }
    
    if (!is.null(kappa3)) {
      if (kappa3 <= 0) {
        stop("Argument 'kappa3' must not be negative.")
      } 
      if (kappa3 > 1) {
        stop("Argument 'kappa3' must not be larger than 1.")
      }
    }
    
    if (kappa4 <= 0) {
      stop("Argument 'kappa4' must not be negative.")
    }
    if (kappa4 > 1) {
      stop("Argument 'kappa4' must not be larger than 1.")
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
          incl_matrix[, (i - 1) * k + 1:k] <- kappa2 / i
          if (k > 1) {
            diag(incl_matrix[, (i - 1) * k + 1:k]) <- kappa1 / i 
          } else {
            incl_matrix[, (i - 1) * k + 1] <- kappa1 / i
          }
        }
      }
      
      if (m > 0) {
        incl_matrix[, n_a + 1:m] <- kappa3
        if (s > 0) {
          for (i in 1:s) {
            incl_matrix[, n_a + m + (i - 1) * m + 1:m] <- kappa3 / (1 + i)
          }
        }
      }
      
      if (n_c > 0) {
        incl_matrix[, n_a + n_b + 1:n_c] <- kappa4
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