
# The function is used within add_initial_values-methods to add initial values
# for the errors of the measurement equation
.add_initial_values_measurement_errors <- function(object, method, u) {
  
  if (method %in% c("ols", "maxlik")) {
    
    y <- t(object[["data"]][["y"]])
    z <- object[["data"]][["z"]]
    k <- nrow(y)
    tt <- ncol(y)
    
    # Errors
    if (object$model$error %in% c("gamma", "gamma+covar")) {
      object[["initial"]][["sigma_i"]] <- diag(1 / apply(u, 1, stats::var), k)
    }
    
    if (object$model$error %in% c("sv", "sv+covar")) {
      u <- apply(u, 1, stats::var)
      object[["initial"]][["h"]] <- log(matrix(u, nrow = NCOL(y), ncol = NROW(y), byrow = TRUE))
      object[["initial"]][["h_init"]] <- matrix(object[["initial"]][["h"]][1, ])
      object[["initial"]][["h_state_variance"]] <- object$priors$sigma$state_variance
      object[["initial"]][["h_offset"]] <- object$priors$sigma$offset
    } 
    
    if (object$model$error == "wishart") {
      object[["initial"]][["sigma_i"]] <- solve(tcrossprod(u) / tt)
    }
  }
  
  if (method == "prior") {
    # Errors
    if (object$model$error %in% c("gamma", "gamma+covar")) {
      sigma_shape <- object$priors$sigma$shape
      sigma_rate <- 1 / object$priors$sigma$rate
      object[["initial"]][["sigma_i"]] <- diag(1, k)
      for (i in 1:k) {
        object[["initial"]][["sigma_i"]][i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
      }
    }
    
    if (object$model$error %in% c("sv", "sv+covar")) {
      mu <- object$priors$sigma$mu
      vinv <- object$priors$sigma$v_i
      h_draw <- mu + chol(vinv) %*% stats::rnorm(NROW(y))
      object[["initial"]][["h"]] <- matrix(h_draw, nrow = NCOL(y), ncol = NROW(y), byrow = TRUE)
      object[["initial"]][["h_init"]] <- matrix(object[["initial"]][["h"]][1, ])
      object[["initial"]][["h_state_variance"]] <- object$priors$sigma$state_variance
      object[["initial"]][["h_offset"]] <- object$priors$sigma$offset
    }
    
    if (object$model$error == "wishart") {
      sigma_df <- object$priors$sigma$df
      sigma_scale <- solve(object$priors$sigma$scale)
      object[["initial"]][["sigma_i"]] <- matrix(stats::rWishart(1, df = sigma_df, Sigma = sigma_scale)[,,1], k)
    }
  }
  
  return(object)
}


# Generates the initial values of the state equation in case of TVP models
.add_initial_values_state_errors <- function(object) {
  
  if (object[["model"]][["tvp"]] & !is.null(object[["data"]][["z"]])) {
    object[["initial"]][["a_v_i"]] <- Matrix::Diagonal(n = ncol(object[["data"]][["z"]]), x = 0)
    for (i in 1:ncol(object[["data"]][["z"]])) {
      object[["initial"]][["a_v_i"]][i, i] <-  1 / stats::rgamma(1,
                                                                 shape = object[["priors"]][["a"]][["shape"]][i],
                                                                 rate = object[["priors"]][["a"]][["rate"]][i])
    }
  }
  
  if (object[["model"]][["tvp"]] & object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & ncol(object[["data"]][["y"]]) > 1) {
    n_psi <- length(object[["priors"]][["psi"]][["shape"]])
    object[["initial"]][["psi_v_i"]] <- Matrix::Diagonal(n = n_psi, x = 0)
    for (i in 1:n_psi) {
      object[["initial"]][["psi_v_i"]][i, i] <-  1 / stats::rgamma(1,
                                                                   shape = object[["priors"]][["psi"]][["shape"]][i],
                                                                   rate = object[["priors"]][["psi"]][["rate"]][i])
    }
  }
  
  return(object)
}
