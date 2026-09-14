#' Stochastic Search Variable Selection Prior
#' 
#' Calculates the priors for a Bayesian VAR model, which employs stochastic search variable selection (SSVS).
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{create_bvarmodel}}
#' or \code{\link{create_bvecmodel}}.
#' @param tau a numeric vector of two elements containing the prior standard errors of restricted
#' variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#' as its second. Default is \code{c(0.05, 10)}.
#' @param semiautomatic an optional numeric vector of two elements containing the factors by which
#' the standard errors associated with an unconstrained least squares estimate of the VAR model are
#' multiplied to obtain the prior standard errors of restricted (\eqn{\tau_0}) and unrestricted
#' (\eqn{\tau_1}) variables. This is the semiautomatic approach described in George et al. (2008).
#' @param ... arguments passed forward to method.
#' 
#' @return A list containing the vectors of prior standard deviations for restricted
#' and unrestricted variables, respectively.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' @examples
#' 
#' # Prepare data
#' data("e6")
#' 
#' # Generate model input
#' object <- create_bvecmodel(e6, p = 2, r = 1)
#' 
#' # Obtain SSVS prior
#' prior <- ssvs_prior(object, semiautomatic = c(.1, 10))
#' 
#' @export
ssvs_prior.bvecmodel <- function(object, tau = c(0.05, 10), semiautomatic = NULL, ...) {
  
  if (object[["model"]][["error"]] %in% c("sv", "sv-covar")) {
    stop("SSVS cannot be used with models with stochastic volatility.")
  }
  
  if (object[["model"]][["tvp"]]) {
    stop("SSVS cannot be used with models with time varying parameter.")
  }
  
  if (!is.null(semiautomatic)) {
    if (!"numeric" %in% class(semiautomatic)) {
      stop("Argument 'semiautomatic' must be a numeric vector of length 2.")
    }
    if (length(semiautomatic) != 2) {
      stop("Argument 'semiautomatic' must be a numeric vector of length 2.")
    }
  }
  
  y <- t(object[["data"]][["train"]][["y"]])
  w <- t(object[["data"]][["train"]][["w"]])
  tt <- NCOL(y)
  k <- NROW(y)

  # The output must have the same length as the number of columns in Z

  if (!is.null(object[["data"]][["train"]][["z"]])) {
    
    tau0 <- NULL
    tau1 <- NULL
    
    # alpha coefficients ----
    if (object$model$rank > 0) {
      tau0 <- rbind(tau0, matrix(1, k * object$model$rank))
      tau1 <- rbind(tau1, matrix(1, k * object$model$rank))
    }
    
    # Non-alpha coefficients ----
    structural <- object[["model"]][["structural"]] & k > 1
    n_struct <- if (structural) k * (k - 1) / 2 else 0
    has_x <- !is.null(object[["data"]][["train"]][["x"]])
    x <- if (has_x) t(object[["data"]][["train"]][["x"]]) else matrix(NA_real_, 0, tt)
    n_x <- nrow(x)

    if (!is.null(semiautomatic) & (has_x | structural)) {

      # As in the VAR method: the least squares standard errors this scales
      # by need more observations than there are regressors per equation, of
      # which the last equation of a structural model has k - 1 more.
      n_max <- n_x + if (structural) k - 1 else 0
      if (tt <= n_max) {
        stop("Argument 'semiautomatic' scales the prior by least squares ",
             "standard errors, but the training sample has ", tt,
             " observations for ", n_max, " regressors per equation. ",
             "Omit 'semiautomatic' to use the fixed values in 'tau', reduce ",
             "the lag order, or provide a longer training sample.")
      }

      if (structural) {
        # Equation i of a structural model regresses the difference of variable
        # i on the regressors and on minus the current differences of the
        # variables before it. The residuals of such a recursive system are
        # orthogonal across equations, so the standard errors are those of least
        # squares equation by equation, each with its own degrees of freedom.
        # The contemporaneous coefficients used to fall back to 'tau'.
        se_x <- matrix(NA_real_, k, n_x)
        se_a0 <- matrix(NA_real_, k, k)
        for (i in 1:k) {
          z_i <- cbind(t(x), -t(y[seq_len(i - 1), , drop = FALSE]))
          if (ncol(z_i) == 0) {
            next
          }
          zz_inv <- solve(crossprod(z_i))
          u_i <- y[i, ] - z_i %*% (zz_inv %*% crossprod(z_i, y[i, ]))
          se_i <- sqrt(diag(zz_inv) * sum(u_i^2) / (tt - ncol(z_i)))
          se_x[i, ] <- se_i[seq_len(n_x)]
          se_a0[i, seq_len(i - 1)] <- se_i[n_x + seq_len(i - 1)]
        }
        se_ols <- c(c(se_x), se_a0[lower.tri(se_a0)])
      } else {
        ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
        u <- y - ols %*% x
        sigma_ols <- tcrossprod(u) / (tt - n_x) # OLS error covariance matrix
        cov_ols <- kronecker(solve(tcrossprod(x)), sigma_ols) # Sqrt of diagonal elements are the t-ratios
        se_ols <- sqrt(diag(cov_ols)) # OLS standard errors
      }

      tau0 <- append(tau0, se_ols * semiautomatic[1]) # Prior if excluded
      tau1 <- append(tau1, se_ols * semiautomatic[2]) # Prior if included

    } else {
      tau0 <- append(tau0, rep(tau[1], k * n_x + n_struct))
      tau1 <- append(tau1, rep(tau[2], k * n_x + n_struct))
    }
    
    result <- list("tau0" = matrix(tau0),
                   "tau1" = matrix(tau1))
    
  } else {
    result <- NULL
  }
  
  return(result)
}