#' Check Gibbs Sampler Input
#' 
#' Checks if the input to function \code{\link{bvecpost}} is suitable for the
#' package's own Gibbs sampling algorithms.
#' 
#' @param object a list containing the necesseary data objects required for
#' a Gibbs sampling algorithm. Usually, the output of a call to 
#' \code{\link{create_vec_model}} in combination with \code{\link{add_priors}}
#' and \code{\link{add_initial_values}}.
#' 
#' 
#' @export
check_bvecpost_input <- function(object){
  
  # Coefficients ----
  if (!is.null(object[["data"]][["z"]])) {
    
    # Priors
    if (is.null(object[["priors"]][["a"]])) {
      stop("Missing element 'object$priors$a'.")
    }
    for (i in c("mu", "v_i")) {
      if (is.null(object[["priors"]][["a"]][[i]])) {
        stop(paste0("Missing element 'object$priors$a$", i, "'."))
      }
    } 
    if (object[["model"]][["tvp"]]) {
      for (i in c("shape", "rate")) {
        if (is.null(object[["priors"]][["a"]][[i]])) {
          stop(paste0("Missing element 'object$priors$a$", i, "'."))
        }
      } 
    }
    
    # Initial values
    if (is.null(object[["initial"]][["a"]])) {
      stop("Missing element 'object$inital$a'.")
    }
    if (object[["model"]][["tvp"]]) {
      for (i in c("a_init", "a_v_i")) {
        if (is.null(object[["initial"]][[i]])) {
          stop(paste0("Missing element 'object$initial$", i, "'."))
        }
      } 
    }
    
  }
  
  # Covariances ----
  if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")) {
    # Priors
    if (is.null(object[["priors"]][["psi"]])) {
      stop("Missing element 'object$priors$psi'.")
    }
    for (i in c("mu", "v_i")) {
      if (is.null(object[["priors"]][["psi"]][[i]])) {
        stop(paste0("Missing element 'object$priors$psi$", i, "'."))
      }
    } 
    if (object[["model"]][["tvp"]]) {
      for (i in c("shape", "rate")) {
        if (is.null(object[["priors"]][["psi"]][[i]])) {
          stop(paste0("Missing element 'object$priors$psi$", i, "'."))
        }
      } 
    }
    
    # Initial values
    if (is.null(object[["initial"]][["psi"]])) {
      stop("Missing element 'object$inital$psi'.")
    }
    if (object[["model"]][["tvp"]]) {
      for (i in c("psi_init", "psi_v_i")) {
        if (is.null(object[["initial"]][[i]])) {
          stop(paste0("Missing element 'object$initial$", i, "'."))
        }
      } 
    }
  }
  
  
  # Errors ----
  if (is.null(object[["priors"]][["sigma"]])) {
    stop("Missing element 'object$priors$sigma'.")
  }
  
  # Wishart prior
  if (object[["model"]][["error"]] == "wishart") {
    # Priors
    for (i in c("df", "scale")) {
      if (is.null(object[["priors"]][["sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    }
  }
  
  # Gamma prior
  if (object[["model"]][["error"]] %in% c("gamma", "gamma+covar")) {
    # Priors
    for (i in c("shape", "rate")) {
      if (is.null(object[["priors"]][["sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    } 
  }
  
  # Stochastic volatility
  if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
    for (i in c("mu", "v_i", "shape", "rate")) {
      if (is.null(object[["priors"]][["sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    } 
  }
}