
# The function is used within add_initial_values-methods to add initial values
# for the errors of the measurement equation
.add_initial_values_measurement_errors <- function(object, method, u) {
  
  if (!method %in% c("ols", "maxlik", "prior")) {
    stop("Unknown specification of argument 'method'.")
  }
  
  if (method %in% c("ols", "maxlik")) {
    
    y <- matrix(t(object[["data"]][["train"]][["y"]]))
    z <- object[["data"]][["train"]][["z"]]
    k <- object[["model"]][["k"]]
    tvp <- object[["model"]][["tvp"]]
    tt <- nrow(y)
    
    # Errors
    if (object[["model"]][["error"]] %in% c("gamma", "gamma+covar")) {
      object[["initial"]][["u_omega_inv"]] <- diag(1 / apply(u, 1, stats::var), k)
    }
    
    if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
      u <- apply(u, 1, stats::var)
      object[["initial"]][["h"]] <- log(matrix(u, nrow = tt, ncol = k, byrow = TRUE))
      object[["initial"]][["h_init"]] <- matrix(object[["initial"]][["h"]][1, ])
      object[["initial"]][["h_state_variance"]] <- object[["priors"]][["sigma"]][["state_variance"]]
      object[["initial"]][["h_offset"]] <- object[["priors"]][["sigma"]][["offset"]]
    } 
    
    if (object[["model"]][["error"]] == "wishart") {
      object[["initial"]][["u_sigma_inv"]] <- solve(tcrossprod(u) / tt)
    }
  }
  
  if (method == "prior") {
    # Errors
    if (object[["model"]][["error"]] %in% c("gamma", "gamma+covar")) {
      sigma_shape <- object[["priors"]][["sigma"]][["shape"]]
      sigma_rate <- 1 / object[["priors"]][["sigma"]][["rate"]]
      object[["initial"]][["u_omega_inv"]] <- diag(1, k)
      for (i in 1:k) {
        object[["initial"]][["u_omega_inv"]][i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i] / 2, rate = sigma_rate[i] / 2)
      }
    }
    
    if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
      mu <- object[["priors"]][["sigma"]][["mu"]]
      vinv <- object[["priors"]][["sigma"]][["v_inv"]]
      h_draw <- mu + chol(vinv) %*% stats::rnorm(tt)
      object[["initial"]][["h"]] <- matrix(h_draw, nrow = tt, ncol = k, byrow = TRUE)
      object[["initial"]][["h_init"]] <- matrix(object[["initial"]][["h"]][1, ])
      object[["initial"]][["h_state_variance"]] <- object[["priors"]][["sigma"]][["state_variance"]]
      object[["initial"]][["h_offset"]] <- object[["priors"]][["sigma"]][["offset"]]
    }
    
    if (object[["model"]][["error"]] == "wishart") {
      sigma_df <- object[["priors"]][["sigma"]][["df"]]
      sigma_scale <- solve(object[["priors"]][["sigma"]][["scale"]])
      object[["initial"]][["u_sigma_inv"]] <- matrix(stats::rWishart(1, df = sigma_df, Sigma = sigma_scale)[,,1], k)
    }
  }
  
  return(object)
}


# Generates the initial values of the state equation in case of TVP models
.add_initial_values_state_errors <- function(object) {
  
  if (object[["model"]][["tvp"]] & !is.null(object[["data"]][["train"]][["z"]])) {
    object[["initial"]][["a_sigma_inv"]] <- diag(0, ncol(object[["data"]][["train"]][["z"]]))
    for (i in 1:ncol(object[["data"]][["train"]][["z"]])) {
      object[["initial"]][["a_sigma_inv"]][i, i] <-  stats::rgamma(1,
                                                                   shape = object[["priors"]][["a"]][["shape"]][i] / 2,
                                                                   rate = object[["priors"]][["a"]][["rate"]][i] / 2)
    }
  }
  
  use_covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  if (object[["model"]][["tvp"]] & use_covar & object[["model"]][["k"]] > 1) {
    n_psi <- length(object[["priors"]][["psi"]][["shape"]])
    object[["initial"]][["psi_sigma_inv"]] <- diag(0, n_psi)
    for (i in 1:n_psi) {
      object[["initial"]][["psi_sigma_inv"]][i, i] <- stats::rgamma(1,
                                                                    shape = object[["priors"]][["psi"]][["shape"]][i] / 2,
                                                                    rate = object[["priors"]][["psi"]][["rate"]][i] / 2)
    }
  }
  
  return(object)
}
