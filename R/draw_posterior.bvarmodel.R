#' Posterior Simulation for Vector Autoregressive Models
#' 
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#' 
#' @param object a list of models, usually, the output of a call to
#' \code{\link{create_var_model}} in combination with \code{\link{add_priors}}
#' and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to each model in argument \code{object}.
#' If \code{NULL} (default), the internal functions \code{\link{bvarpost}} is used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of the class of the output of the applied posterior
#' simulation function. In case the package's own function is used, this will
#' result in an object of class \code{"bvar"}.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_var_model(e1, p = 2, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior.bvarmodel <- function(object, posterior_function = NULL, ...){
  
  model <- object
  if ("posteriors" %in% names(model)) {
    model[["posteriors"]] <- NULL
  }
  
  if (is.null(posterior_function)) {
    object <- try(bvarpost(object))
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  } else {
    
    if (is.null(posterior_function)) {
      
      k <- object[["model"]][["k"]]
      p <- object[["model"]][["p"]]
      m <- object[["model"]][["m"]]
      s <- object[["model"]][["s"]]
      n_a <- k * k * p
      n_b <- k * m * s
      n_c <- k * object[["model"]][["n"]]
      tt <- NROW(object[["data"]][["y"]])
      n_z <- ncol(object[["data"]][["z"]])
      draws <- object[["model"]][["iterations"]]
      tvp <- object[["model"]][["tvp"]]
      
      # Modelled variables ----
      draws_a <- NULL
      if (p > 0) {
        pos_a <- 1:n_a
        if (tvp) {
          pos_a <- rep(pos_a, tt) + rep(0:(tt - 1), each = length(pos_a)) * n_z
        }
        draws_a[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_a])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          draws_a[["sigma"]] <- matrix(NA, n_a, draws)
          draws_a[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_a[1:n_a]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          draws_a[["lambda"]] <- matrix(NA, n_a, draws)
          draws_a[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_a[1:n_a]])
        }
      }
      
      # Unmodelled variables ----
      draws_b <- NULL
      if (m > 0) {
        pos_b <- n_a + 1:n_b
        if (tvp) {
          pos_b <- rep(pos_b, tt) + rep(0:(tt - 1), each = length(pos_b)) * n_z
        }
        draws_b[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_b])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          draws_b[["sigma"]] <- matrix(NA, n_b, draws)
          draws_b[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_b[1:n_b]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          draws_b[["lambda"]] <- matrix(NA, n_b, draws)
          draws_b[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_b[1:n_b]])
        }
      }
      
      # Deterministic terms ----
      draws_c <- NULL
      if (n_c > 0) {
        pos_c <- n_a + n_b + 1:n_c
        if (tvp) {
          pos_c <- rep(pos_c, tt) + rep(0:(tt - 1), each = length(pos_c)) * n_z
        }
        draws_c[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_c])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          draws_c[["sigma"]] <- matrix(NA, n_c, draws)
          draws_c[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_c[1:n_c]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          draws_c[["lambda"]] <- matrix(NA, n_c, draws)
          draws_c[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_c[1:n_c]])
        }
      }
      
      # Structural ----
      draws_a0 <- .extract_posterior_a0(object)
      
      # Sigma----
      draws_sigma <- .extract_posterior_sigma(object)
      
      # Create bvar object
      object <- bvar(data = NULL,
                     exogen = NULL,
                     y = object[["data"]][["y"]],
                     x = object[["data"]][["x"]],
                     A0 = draws_a0,
                     A = draws_a,
                     B = draws_b,
                     C = draws_c,
                     Sigma = draws_sigma)
    }
  }
  
  return(object)
}