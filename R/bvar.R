#' Bayesian Vector Autoregression Objects
#' 
#' \code{bvar} is used to create objects of class \code{"bvar"}.
#' 
#' @param data the original time-series object of endogenous variables.
#' @param exogen the original time-series object of unmodelled variables.
#' @param y a time-series object of endogenous variables,
#' usually, a result of a call to \code{\link{gen_var}}.
#' @param x a time-series object of \eqn{(pK + (1+s)M + N)} regressor variables, usually, a result of a
#' call to \code{\link{gen_var}}.
#' @param A0 either a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC coefficient
#' draws of structural parameters and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param A either a \eqn{pK^2 \times S} matrix of MCMC coefficient draws of lagged endogenous variables or
#' a named list, where element \code{coeffs} contains a \eqn{pK^2 \times S} matrix of MCMC coefficient draws
#' of lagged endogenous variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param B either a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables
#' or a named list, where element \code{coeffs} contains a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of
#' unmodelled, non-deterministic variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param C either a \eqn{KN \times S} matrix of MCMC coefficient draws of deterministic terms or
#' a named list, where element \code{coeffs} contains a \eqn{KN \times S} matrix of MCMC coefficient draws of
#' deterministic terms and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param Sigma a \eqn{K^2 \times S} matrix of MCMC draws for the error variance-covariance matrix or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC draws for the
#' error variance-covariance matrix and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed to the covariances.
#' 
#' @details For the VARX model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_i y_{t-i} + \sum_{i = 0}^{s} B_i x_{t - i} + C d_t + u_t}
#' the function collects the S draws of a Gibbs sampler (after the burn-in phase) in a standardised object,
#' where \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of lagged endogenous variabels.
#' \eqn{x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{B_i} its corresponding coefficient matrix.
#' \eqn{d_t} is an N-dimensional vector of deterministic terms
#' and \eqn{C} its corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#' 
#' The draws of the different coefficient matrices provided in \code{A0}, \code{A},
#' \code{B}, \code{C} and \code{Sigma} have to correspond to the same MCMC iterations.
#' 
#' @return An object of class \code{"bvar"} containing the following components, if specified:
#' \item{data}{the original time-series object of endogenous variables.}
#' \item{exogen}{the original time-series object of unmodelled variables.}
#' \item{y}{a \eqn{K \times T} matrix of endogenous variables.}
#' \item{x}{a \eqn{(pK + (1+s)M + N) \times T} matrix of regressor variables.}
#' \item{A0}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of structural parameters.}
#' \item{A0_lambda}{an \eqn{S \times K^2} "mcmc" object of inclusion parameters for structural parameters.}
#' \item{A}{an \eqn{S \times pK^2} "mcmc" object of coefficient draws of lagged endogenous variables.}
#' \item{A_lambda}{an \eqn{S \times pK^2} "mcmc" object of inclusion parameters for lagged endogenous variables.}
#' \item{B}{an \eqn{S \times ((1 + s)MK)} "mcmc" object of coefficient draws of unmodelled, non-deterministic variables.}
#' \item{B_lambda}{an \eqn{S \times ((1 + s)MK)} "mcmc" object of inclusion parameters for unmodelled, non-deterministic variables.}
#' \item{C}{an \eqn{S \times NK} "mcmc" object of coefficient draws of deterministic terms.}
#' \item{C_lambda}{an \eqn{S \times NK} "mcmc" object of inclusion parameters for deterministic terms.}
#' \item{Sigma}{an \eqn{S \times K^2} "mcmc" object of variance-covariance draws.}
#' \item{Sigma_lambda}{an \eqn{S \times K^2} "mcmc" object of inclusion parameters for error covariances.}
#' \item{specifications}{a list containing information on the model specification.}

#' @examples
#' 
#' # Get data
#' data("e1")
#' e1 <- diff(log(e1))
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' data <- gen_var(e1, p = 2, deterministic = "const")
#'
#' # Add priors
#' model <- add_priors(data,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 0, scale = .00001))
#' 
#' # Set RNG seed for reproducibility 
#' set.seed(1234567)
#' 
#' iterations <- 400 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations and burnin should be much higher.
#' burnin <- 100 # Number of burn-in draws
#' draws <- iterations + burnin # Total number of MCMC draws
#'
#' y <- t(model$data$Y)
#' x <- t(model$data$Z)
#' tt <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' m <- k * nrow(x) # Number of estimated coefficients
#' 
#' # Priors
#' a_mu_prior <- model$priors$coefficients$mu # Vector of prior parameter means
#' a_v_i_prior <- model$priors$coefficients$v_i # Inverse of the prior covariance matrix
#' 
#' u_sigma_df_prior <- model$priors$sigma$df # Prior degrees of freedom
#' u_sigma_scale_prior <- model$priors$sigma$scale # Prior covariance matrix
#' u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
#'
#' # Initial values
#' u_sigma_i <- diag(1 / .00001, k)
#'
#' # Data containers for posterior draws
#' draws_a <- matrix(NA, m, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#'
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'  # Draw conditional mean parameters
#'  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
#'
#'  # Draw variance-covariance matrix
#'  u <- y - matrix(a, k) %*% x # Obtain residuals
#'  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
#'  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'
#'  # Store draws
#'  if (draw > burnin) {
#'   draws_a[, draw - burnin] <- a
#'   draws_sigma[, draw - burnin] <- solve(u_sigma_i)
#'  }
#' }
#' 
#' # Generate bvar object
#' bvar_est <- bvar(y = model$data$Y, x = model$data$Z,
#'                  A = draws_a[1:18,], C = draws_a[19:21, ],
#'                  Sigma = draws_sigma)
#'                  
#' @export
bvar <- function(data = NULL, exogen = NULL, y = NULL, x = NULL,
                 A0 = NULL, A = NULL, B = NULL,
                 C = NULL, Sigma = NULL) {

  result <- NULL
  if (is.null(y)) {
    stop("At least argument 'y' must be specified.")
  } else {
    result$y <- y
    k <- NCOL(y)
  }
  
  if(!is.null(A)) {
    if ("list" %in% class(A)) {
      result$A <- coda::mcmc(t(A$coeffs))
      result$A_lambda <- coda::mcmc(t(A$lambda))
    } else {
      result$A <- coda::mcmc(t(A)) 
    }
    n_a <- ncol(result$A)
    if (n_a %% k == 0) {
      p <- n_a / k^2 
    } else {
      stop("Row number of argument 'A' is not a multiple of the number of endogenous variables.")
    }
  } else {
    p <- 0
  }
  
  if(!is.null(data)) {
    result$data <- data
  }
  if(!is.null(exogen)) {
    result$exogen <- exogen
  }
  if(!is.null(x)) {
    result$x <- x
  }
  if(!is.null(A0)) {
    if ("list" %in% class(A0)) {
      result$A0 <- coda::mcmc(t(A0$coeffs))
      result$A0_lambda <- coda::mcmc(t(A0$lambda))
    } else {
      result$A0 <- coda::mcmc(t(A0)) 
    }
  }
  if(!is.null(B)) {
    if ("list" %in% class(B)) {
      result$B <- coda::mcmc(t(B$coeffs))
      result$B_lambda <- coda::mcmc(t(B$lambda))
    } else {
      result$B <- coda::mcmc(t(B)) 
    }
  }
  if(!is.null(C)) {
    if ("list" %in% class(C)) {
      result$C <- coda::mcmc(t(C$coeffs))
      result$C_lambda <- coda::mcmc(t(C$lambda))
    } else {
      result$C <- coda::mcmc(t(C)) 
    }
  }
  if(!is.null(Sigma)) {
    if ("list" %in% class(Sigma)) {
      result$Sigma <- coda::mcmc(t(Sigma$coeffs))
      result$Sigma_lambda <- coda::mcmc(t(Sigma$lambda))
    } else {
      result$Sigma <- coda::mcmc(t(Sigma)) 
    }
  }
  
  result$specifications <- list("dims" = c("K" = k),
                                "lags" = c("p" = p))
  
  if(!is.null(B)) {
    if (!is.null(exogen)) {
      m <- NCOL(exogen)
      n_b <- ncol(result$B)
      if (n_b %% k == 0) {
        s <- n_b / (m * k) 
      } else {
        stop("Row number of argument 'B' is not a multiple of the number of endogenous variables.")
      }
      result$specifications$dims["M"] <- m
      result$specifications$lags["s"] <- s
    }
  }
  
  class(result) <- "bvar"
  return(result)
}