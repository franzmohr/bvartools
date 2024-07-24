#' Thinning Posterior Draws
#'
#' Thins the MCMC posterior draws in an object of class \code{"dfm"}.
#'
#' @param x an object of class \code{"dfm"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
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
#' object <- draw_posterior(model)
#'
#' # Plot factors
#' object <- thin(object, thin = 2)
#'
#' @return An object of class \code{"dfm"}.
#'
#' @export
thin.dfm <- function(x, thin = 10, ...) {
  
  # Detect number of posterior draws
  draws <- NA
  if (!is.null(x[["lambda"]])) {
    draws <- nrow(x[["lambda"]])
  }
  vars <- c("u", "v", "a")
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[[i]])) {
        draws <- nrow(x[[i]])
      }
    }
  }
  
  # Determine the kept observations
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  # Thinning
  vars <- c("lambda", "factor", "a", "u", "v")
  for (i in vars) {
    if (!is.null(x[[i]])) {
      x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin)
    }
  }
  
  return(x)
}
