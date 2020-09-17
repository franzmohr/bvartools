#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvar"} or \code{"bvec"}.
#' 
#' @param x an object of class \code{"bvar"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Obtain data matrices
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' object <- thin_posterior(object)
#' 
#' @return An object of class \code{"bvar"}.
#' 
#' @export
thin_posterior.bvar <- function(x, thin = 10) {
  
  draws <- NA
  if (!is.null(x[["A"]])) {
    draws <- nrow(x[["A"]])
  }
  vars <- c("B", "C", "A0")
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
  
  vars <- c("A", "B", "C", "A0", "Sigma")
  
  for (i in vars) {
    if (!is.null(x[[i]])) {
      x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin) 
    }
  }
  
  return(x)
}