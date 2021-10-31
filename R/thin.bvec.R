#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvec"}.
#' 
#' @param x an object of class \code{"bvec"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples 
#' 
#' # Load data
#' data("e6")
#' 
#' # Generate model data
#' model <- gen_vec(e6, p = 2, r = 1,
#'                  const = "unrestricted", seasonal = "unrestricted",
#'                  iterations = 100, burnin = 10)
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Thin
#' object <- thin(object)
#' 
#' @return An object of class \code{"bvec"}.
#' 
#' @export
thin.bvec <- function(x, thin = 10, ...) {
  
  draws <- NA
  vars <- c("Pi", "Pi_x", "Pi_d", "Gamma", "Upsilon", "C", "A0")
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[[i]])) {
        if (x[["specifications"]][["tvp"]][[i]]) {
          draws <- nrow(x[[i]][[1]])
        } else {
          draws <- nrow(x[[i]]) 
        }
      }   
    }
  }
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  vars <- c("alpha",
            "beta", "beta_x", "beta_d",
            "Pi", "Pi_x", "Pi_d",
            "Gamma", "Gamma_sigma", "Gamma_lambda",
            "Upsilon", "Upsilon_sigma", "Upsilon_lambda",
            "C", "C_sigma", "C_lambda",
            "A0", "A0_sigma", "A0_lambda",
            "Sigma")
  for (i in vars) {
    if (!is.null(x[[i]])) {
      if (is.list(x[[i]])) {
        for (j in 1:length(x[[i]])) {
          x[[i]][[j]] <- coda::mcmc(as.matrix(x[[i]][[j]][pos_thin,]), start = start, end = end, thin = thin)  
        }
      } else {
        x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin)  
      }
    }
  }
  
  return(x)
}