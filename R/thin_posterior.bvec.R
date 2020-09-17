#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvec"}.
#' 
#' @param x an object of class \code{"bvec"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
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
#' object <- thin_posterior(object)
#' 
#' @return An object of class \code{"bvec"}.
#' 
#' @export
thin_posterior.bvec <- function(x, thin = 10) {
  
  draws <- NA
  if (!is.null(x[["Pi"]])) {
    draws <- nrow(x[["Pi"]])
  }
  vars <- c("Pi_x", "Pi_d", "Gamma", "Upsilon", "C", "A0")
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[[i]])) {
        draws <- nrow(x[[i]])
      }
    }
  }
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  vars <- c("alpha", "beta", "Pi", "Pi_x", "Pi_d", "Gamma", "Upsilon", "C", "A0", "Sigma")
  for (i in vars) {
    if (!is.null(x[[i]])) {
      x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin) 
    }
  }
  
  return(x)
}