#' Check Gibbs Sampler Input
#' 
# Checks if the input to function \code{\link{bvarpost}} is suitable for the
# package's own Gibbs sampling algorithms.
# 
# @param object a list containing the necesseary data objects required for
# a Gibbs sampling algorithm. Usually, the output of a call to 
# \code{\link{create_bvarmodel}} in combination with \code{\link{add_priors}}
# and \code{\link{add_initial_values}}.
# 
# 
# @export
.check_bvarpost_input <- function(object){
  
  # Coefficients ----
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    
    # Priors
    if (is.null(object[["priors"]][["a"]])) {
      stop("Missing element 'object$priors$a'.")
    }
    for (i in c("mu", "v_inv")) {
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
      for (i in c("a_init", "a_sigma_inv")) {
        if (is.null(object[["initial"]][[i]])) {
          stop(paste0("Missing element 'object$initial$", i, "'."))
        }
      } 
    }
    
  }
  
  # Covariances ----
  if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & object[["model"]][["k"]] > 1) {
    # Priors
    if (is.null(object[["priors"]][["psi"]])) {
      stop("Missing element 'object$priors$psi'.")
    }
    for (i in c("mu", "v_inv")) {
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
      for (i in c("psi_init", "psi_sigma_inv")) {
        if (is.null(object[["initial"]][[i]])) {
          stop(paste0("Missing element 'object$initial$", i, "'."))
        }
      } 
    }
  }
  
  
  # Errors ----
  if (is.null(object[["priors"]][["u_sigma"]])) {
    stop("Missing element 'object$priors$u_sigma'.")
  }
  
  # Wishart prior
  if (object[["model"]][["error"]] == "wishart") {
    # Priors
    for (i in c("df", "scale")) {
      if (is.null(object[["priors"]][["u_sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    }
  }
  
  # Gamma prior
  if (object[["model"]][["error"]] %in% c("gamma", "gamma+covar")) {
    # Priors
    for (i in c("shape", "rate")) {
      if (is.null(object[["priors"]][["u_sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    } 
  }
  
  # Stochastic volatility
  if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
    for (i in c("mu", "v_inv", "shape", "rate")) {
      if (is.null(object[["priors"]][["u_sigma"]][[i]])) {
        stop(paste0("Missing element 'object$priors$sigma$", i, "'."))
      }
    } 
  }
}