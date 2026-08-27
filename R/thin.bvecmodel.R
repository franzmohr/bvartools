#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class 'bvecmodel'.
#' 
#' @param x an object of class 'bvecmodel'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_bvecmodel(e6, p = 2, r = 1,
#'                          const = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Thinning
#' model <- thin(model, 2)
#' 
#' @export
#' @method thin bvecmodel
thin.bvecmodel <- function(x, thin = 10, ...) {
  
  draws <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  # Posteriors with sub-lists
  vars <- c("a", "beta", "psi", "u_sigma_inv", "u_omega_inv")
  subvars <- c("coeffs", "lambda")
  for (i in vars) {
    if (!is.null(x[["posterior"]][[i]])){
      for (j in subvars) {
        if (!is.null(x[["posterior"]][[i]][[j]])) {
          x[["posterior"]][[i]][[j]] <- coda::mcmc(as.matrix(x[["posterior"]][[i]][[j]][pos_thin,]), start = start, end = end, thin = thin)
        }
      } 
    }
  }
  
  # Posteriors w/o sub-lists
  vars <- c("loglik", "forecast", "forecast_errors")
  for (i in vars) {
    if (!is.null(x[["posterior"]][[i]])){
      x[["posterior"]][[i]] <- coda::mcmc(as.matrix(x[["posterior"]][[i]][pos_thin,]), start = start, end = end, thin = thin)
    }
  }

  return(x)
}