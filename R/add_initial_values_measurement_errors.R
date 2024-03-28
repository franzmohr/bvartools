
# The function is used within add_initial_values-methods to add initial values
# for the errors of the measurement equation
.add_initial_values_measurement_errors <- function(object, method, u) {
  
  if (method == "ols") {
    
    y <- t(object[["data"]][["Y"]])
    z <- object[["data"]][["SUR"]]
    k <- nrow(y)
    tt <- ncol(y)
    if (!is.null(z)) {
      est_vars <- ncol(z)
    } else {
      est_vars <- 0
    }
    
    # Errors
    if (object[["model"]][["sv"]]) {
      
      stop("Initial values for SV not implemented yet.")
      
      u <- apply(u, 1, stats::var)
      object[["initial"]][["sigma"]][["h"]] <- log(matrix(u, nrow = NCOL(y), ncol = NROW(y), byrow = TRUE))
      object[["initial"]][["sigma"]][["sigma_h"]] <- matrix(sigma[["sigma_h"]], NROW(y))
      
      if (is.null(sigma[["constant"]])) {
        warning("Argument 'sigma$constant' not specified. Using the value 0.0001.")
        sigma[["constant"]] <- .0001
      }
      object[["initial"]][["sigma"]][["constant"]] <- matrix(sigma[["constant"]], NROW(y))
      
    } else {
      
      if (object$priors$V$type == "wishart") {
        object[["initial"]][["V_i"]] <- solve(tcrossprod(u) / tt)
      }
      if (object$priors$V$type == "gamma") {
        object[["initial"]][["V_i"]] <- diag(1 / apply(u, 1, stats::var), k)
      }
      
    }
  }
  
  if (method == "prior") {
    # Errors
    if (object[["model"]][["sv"]]) {
      
      stop("Initial values for SV not implemented yet.")
      
      if (is.null(sigma[["constant"]])) {
        warning("Argument 'sigma$constant' not specified. Using the value 0.0001.")
        sigma[["constant"]] <- .0001
      }
      object[["initial"]][["sigma"]][["constant"]] <- matrix(sigma[["constant"]], NROW(y))
      
    } else {
      
      if (object$priors$V$type == "wishart") {
        sigma_df <- object$priors$V$df
        sigma_scale <- solve(object$priors$V$scale)
        object[["initial"]][["V_i"]] <- matrix(stats::rWishart(1, df = sigma_df, Sigma = sigma_scale)[,,1], k)
      }
      if (object$priors$V$type == "gamma") {
        sigma_shape <- object$priors$V$shape
        sigma_rate <- 1 / object$priors$V$rate
        object[["initial"]][["V_i"]] <- diag(1, k)
        for (i in 1:k) {
          object[["initial"]][["V_i"]][i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
        }
      }
    }
  }
  
  return(object)
}