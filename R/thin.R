#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvar"} or \code{"bvec"}.
#' 
#' @param object an object of class \code{"bvar"} or \code{"bvec"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' 
#' @examples 
#' data("e6")
#' data <- gen_vec(e6, p = 4, const = "unrestricted", season = "unrestricted")
#' 
#' y <- data$Y
#' w <- data$W
#' x <- data$X
#' 
#' # Reset random number generator for reproducibility
#' set.seed(1234567)
#' 
#' iter <- 4000 # Number of iterations of the Gibbs sampler
#' burnin <- 1000 # Number of burn-in draws
#' store <- iter - burnin
#' 
#' r <- 1 # Set rank
#' 
#' t <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' k_w <- nrow(w) # Number of regressors in error correction term
#' k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
#' 
#' k_alpha <- k * r # Number of elements in alpha
#' k_beta <- k_w * r # Number of elements in beta
#' k_gamma <- k * k_x
#' 
#' # Set uninformative priors
#' a_mu_prior <- matrix(0, k_x * k) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, k_x * k) # Inverse of the prior covariance matrix
#' 
#' v_i <- 0
#' p_tau_i <- diag(1, k_w)
#' 
#' u_sigma_df_prior <- r # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
#' u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom
#' 
#' # Initial values
#' beta <- matrix(c(1, -4), k_w, r)
#' 
#' u_sigma_i <- diag(.0001, k)
#' u_sigma <- solve(u_sigma_i)
#' 
#' g_i <- u_sigma_i
#' 
#' # Data containers
#' draws_alpha <- matrix(NA, k_alpha, store)
#' draws_beta <- matrix(NA, k_beta, store)
#' draws_pi <- matrix(NA, k * k_w, store)
#' draws_gamma <- matrix(NA, k_gamma, store)
#' draws_sigma <- matrix(NA, k^2, store)
#' 
#' # Start Gibbs sampler
#' for (draw in 1:iter) {
#'   # Draw conditional mean parameters
#'   temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = u_sigma_i,
#'                          v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
#'                          gamma_mu_prior = a_mu_prior,
#'                          gamma_V_i_prior = a_v_i_prior)
#'   alpha <- temp$alpha
#'   beta <- temp$beta
#'   Pi <- temp$Pi
#'   gamma <- temp$Gamma
#'   
#'   # Draw variance-covariance matrix
#'   u <- y - Pi %*% w - matrix(gamma, k) %*% x
#'   u_sigma_scale_post <- solve(tcrossprod(u) +
#'      v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
#'   u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'   u_sigma <- solve(u_sigma_i)
#'   
#'   # Update g_i
#'   g_i <- u_sigma_i
#'   
#'   # Store draws
#'   if (draw > burnin) {
#'     draws_alpha[, draw - burnin] <- alpha
#'     draws_beta[, draw - burnin] <- beta
#'     draws_pi[, draw - burnin] <- Pi
#'     draws_gamma[, draw - burnin] <- gamma
#'     draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#' 
#' # Number of non-deterministic coefficients
#' k_nondet <- (k_x - 4) * k
#' 
#' # Generate bvec object
#' bvec_est <- bvec(y = y, w = w, x = x,
#'                  Pi = draws_pi,
#'                  Gamma = draws_gamma[1:k_nondet,],
#'                  C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
#'                  Sigma = draws_sigma)
#' 
#' # Thin posterior draws
#' bvec_est <- thin(bvec_est, thin = 5)
#' 
#' @return An object of class \code{"bvar"} or \code{"bvec"}.
#' 
#' @export
thin <- function(object, thin = 5) {
  if (!any(class(object) %in% c("bvar", "bvec"))) {
    stop("Argument 'object' must be of class 'bvar' or 'bvec'.")
  }
  
  if (any(class(object) == "bvar")) {
    var <- TRUE
  }
  if (any(class(object) == "bvec")) {
    var <- FALSE
  }
  
  if (var) {
    draws <- nrow(object$A)
  } else {
    draws <- nrow(object$Pi)
  }
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  if (var) {
    object$A <- coda::mcmc(object$A[pos_thin,], start = start, end = end, thin = thin)
    if(!is.null(object$A0)) {
      object$A0 <- coda::mcmc(object$A0[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$B)) {
      object$B <- coda::mcmc(object$B[pos_thin,], start = start, end = end, thin = thin)
    }
  } else {
    object$Pi <- coda::mcmc(object$Pi[pos_thin,], start = start, end = end, thin = thin)
    if(!is.null(object$A0)) {
      object$A0 <- coda::mcmc(object$A0[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$alpha)) {
      object$alpha <- coda::mcmc(object$alpha[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$beta)) {
      object$beta <- coda::mcmc(object$beta[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Pi_x)) {
      object$Pi_x <- coda::mcmc(object$Pi_x[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Pi_d)) {
      object$Pi_d <- coda::mcmc(object$Pi_d[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Gamma)) {
      object$Gamma <- coda::mcmc(object$Gamma[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Upsilon)) {
      object$Upsilon <- coda::mcmc(object$Upsilon[pos_thin,], start = start, end = end, thin = thin)
    }
  }
  if(!is.null(object$C)) {
    object$C <- coda::mcmc(object$C[pos_thin,], start = start, end = end, thin = thin)
  }
  if(!is.null(object$Sigma)) {
    object$Sigma <- coda::mcmc(object$Sigma[pos_thin,], start = start, end = end, thin = thin)
  }
  return(object)
}