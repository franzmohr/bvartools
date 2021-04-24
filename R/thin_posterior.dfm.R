#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"dfm"}.
#' 
#' @param x an object of class \code{"dfm"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' 
#' @examples 
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(v_i = .01),
#'                     sigma_u = list(shape = 5, rate = 4),
#'                     a = list(v_i = .01),
#'                     sigma_v = list(shape = 5, rate = 4))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Plot factors
#' object <- thin_posterior(object, thin = 2)
#' 
#' @return An object of class \code{"dfm"}.
#' 
#' @export
thin_posterior.dfm <- function(x, thin = 10) {
  
  draws <- NA
  if (!is.null(x[["lambda"]])) {
    draws <- nrow(x[["lambda"]])
  }
  vars <- c("a", "sigma_u", "sigma_v")
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
  
  vars <- c("lambda", "a", "sigma_u", "sigma_v")
  
  for (i in vars) {
    if (!is.null(x[[i]])) {
      x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin) 
    }
  }
  
  return(x)
}