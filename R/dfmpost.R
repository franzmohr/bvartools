#' Posterior Simulation for Dynamic Factor Models
#'
#' Produces draws from the posterior distributions of Bayesian dynamic factor models.
#'
#' @param object an object of class \code{"dfmodel"}, usually, a result of a
#' call to \code{\link{create_df_model}} in combination with \code{\link{add_priors}}
#' and \code{\link{add_initial_values}}.
#'
#' @details The function implements the posterior simulation algorithm for Bayesian dynamic factor models.
#'
#' The implementation follows the description in Chan et al. (2019) and C++ is used to reduce calculation time.
#'
#' @return An object of class \code{"dfm"}.
#'
#' @references
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_df_model(x = bem_dfmdata, p = 1, n = 1,
#'                          iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' object <- dfmpost(model)
#'
#' @export
dfmpost <- function(object) {
  
  # General specifications
  m <- object$model$m
  n_factors <- object$model$n
  p <- object$model$p
  use_a <- p > 0
  if (use_a) {
    n_a <- n_factors * n_factors * p
  }
  tt <- nrow(object$data$x)
  
  # ****************************************************************************
  # Priors ----
  lambda_prior_vinv <- object$priors$lambda$vinv
  
  u_prior_shape <- object$priors$u$shape
  u_prior_rate <- object$priors$u$rate
  
  v_prior_shape <- object$priors$v$shape
  v_prior_rate <- object$priors$v$rate
  
  if (use_a) {
    a_prior_mu <- object$priors$a$mu
    a_prior_vinv <- object$priors$a$vinv
  }
  
  # ****************************************************************************
  # Initial values ----
  
  xvec <- matrix(t(object$data$x))
  
  diag_tt <- Matrix::Diagonal(tt)
  
  lambda <- as.matrix(diag(1, m)[, 1:n_factors])
  lambda[lower.tri(lambda)] <- object$initial$lambda
  
  uinv <- Matrix::Matrix(object$initial$uinv)
  vinv <- object$initial$vinv
  
  H_a <- Matrix::Diagonal(tt * n_factors)
  if (use_a) {
    a <- object$initial$a
    a_mat <- matrix(a, n_factors)
    a_mat_t <- matrix(NA, n_factors * p, n_factors)
    z_a <- matrix(0, tt * n_factors, p * n_factors * n_factors)
  }
  
  iterations <- object$model$iterations
  burnin <- object$model$burnin
  draws <- iterations + burnin
  
  # Storage objects
  draws_lambda <- matrix(NA, m * n_factors, iterations)
  draws_fac <- matrix(NA, tt * n_factors, iterations)
  draws_u <- matrix(NA, m, iterations)
  draws_v <- matrix(NA, n_factors, iterations)
  if (use_a) {
    draws_a <- matrix(NA, n_a, iterations)
  } else {
    draws_a <- NULL
  }
  
  for (draw in 1:draws) {
    
    # Update H_a
    if (use_a) {
      H_a <- bvartools::generate_lower_block_diagonal(a, n_factors, tt)
    }
    
    # Draw factor
    K_f <- crossprod(H_a, kronecker(diag_tt, vinv)) %*% H_a + kronecker(diag_tt, crossprod(lambda, uinv) %*% lambda)
    f_hat <- solve(K_f, kronecker(diag_tt, crossprod(lambda, uinv)) %*% xvec)
    f <- matrix(f_hat + solve(chol(K_f), stats::rnorm(tt * n_factors)))
    ff <- matrix(f, n_factors)
    
    # Draw lambda equation by equation
    lambda <- .post_lambda(object$data$x, ff, prior_vinv = lambda_prior_vinv, uinv = uinv, lambda = lambda)
    
    # Draw uinv
    u = xvec - matrix(lambda %*% ff)
    uinv <- bvartools::post_gamma_measurement_variance(u, u_prior_shape, u_prior_rate, inverse = TRUE)
    
    # Draw vinv
    v <- matrix(H_a %*% f)
    vinv <- bvartools::post_gamma_measurement_variance(v, v_prior_shape, v_prior_rate, inverse = TRUE)
    
    # Draw a
    if (use_a) {
      # Obtain data object with zeros as pre-sample values
      x_a <- matrix(0, n_factors * p, tt)
      for (i in 1:p) {
        x_a[(i - 1) + 1:n_factors, -(1:i)] <- ff[, -((tt - i + 1):tt)]
      }
      z_a <- kronecker(t(x_a), diag(1, n_factors))
      
      K_a <- a_prior_vinv + crossprod(z_a, kronecker(diag_tt, vinv) %*% z_a)
      a_hat <- solve(K_a, a_prior_vinv %*% a_prior_mu + crossprod(z_a, kronecker(diag_tt, vinv) %*% f))
      a <- matrix(a_hat + solve(chol(K_a), stats::rnorm(n_a)))
    }
    
    if (draw > burnin) {
      pos_draw <- draw - burnin
      draws_lambda[, pos_draw] <- matrix(lambda)
      draws_fac[, pos_draw] <- f
      draws_u[, pos_draw] <- 1 / Matrix::diag(uinv)
      draws_v[, pos_draw] <- 1 / Matrix::diag(vinv)
      if (use_a) {
        draws_a[, pos_draw] <- a
      }
    }
  }
  
  object <- dfm(x = object$data$x,
                lambda = draws_lambda,
                fac = draws_fac,
                u = draws_u,
                a = draws_a,
                v = draws_v)
  
  return(object)
}
