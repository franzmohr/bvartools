#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvar"}.
#' 
#' @param x an object of class \code{"bvar"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
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
#' object <- thin(object)
#' 
#' @return An object of class \code{"bvar"}.
#' 
#' @export
thin.bvar <- function(x, thin = 10, ...) {
  
  draws <- NA
  vars <- c("A0", "A", "B", "C", "Sigma")
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
  
  vars <- c("A", "A_sigma", "A_lambda",
            "B", "B_sigma", "B_lambda",
            "C", "C_sigma", "C_lambda",
            "A0", "A0_sigma", "A0_lambda",
            "Sigma", "Sigma_sigma", "Sigma_lambda")
  
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