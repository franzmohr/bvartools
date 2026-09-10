
# The function is used within add_initial_values-methods to add initial values
# for the errors of the measurement equation
.add_initial_values_measurement_errors <- function(object, method, u) {
  
  if (!method %in% c("ols", "maxlik", "prior")) {
    stop("Unknown specification of argument 'method'.")
  }
  
  # The calling method hands over the residuals as a k x T matrix, which is the
  # only place where the number of observations can be read off reliably. The
  # endogenous variables are stored in SUR form, so their number of rows is
  # k * T and not T.
  k <- object[["model"]][["k"]]
  tt <- ncol(u)
  error <- object[["model"]][["error"]]
  sigma_prior <- object[["priors"]][["u_sigma"]]
  
  if (method %in% c("ols", "maxlik")) {
    
    # Errors
    if (error %in% c("gamma", "gamma+covar")) {
      object[["initial"]][["u_omega_inv"]] <- diag(1 / apply(u, 1, stats::var), k)
    }
    
    if (error %in% c("sv", "sv+covar")) {
      h <- log(matrix(apply(u, 1, stats::var), nrow = tt, ncol = k, byrow = TRUE))
      object[["initial"]][["h"]] <- h
      object[["initial"]][["h_init"]] <- matrix(h[1, ])
    } 
    
    if (error == "wishart") {
      object[["initial"]][["u_sigma_inv"]] <- solve(tcrossprod(u) / tt)
    }
    
    if (error == "ald") {
      object <- .add_initial_values_ald(object, u = u, from_prior = FALSE)
    }
  }
  
  if (method == "prior") {
    # Errors
    if (error %in% c("gamma", "gamma+covar")) {
      sigma_shape <- sigma_prior[["shape"]]
      sigma_rate <- 1 / sigma_prior[["rate"]]
      object[["initial"]][["u_omega_inv"]] <- diag(1, k)
      for (i in 1:k) {
        object[["initial"]][["u_omega_inv"]][i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i] / 2, rate = sigma_rate[i] / 2)
      }
    }
    
    if (error %in% c("sv", "sv+covar")) {
      mu <- sigma_prior[["mu"]]
      vinv <- sigma_prior[["v_inv"]]
      h_draw <- mu + chol(vinv) %*% stats::rnorm(k)
      h <- matrix(h_draw, nrow = tt, ncol = k, byrow = TRUE)
      object[["initial"]][["h"]] <- h
      object[["initial"]][["h_init"]] <- matrix(h[1, ])
    }
    
    if (error == "wishart") {
      sigma_df <- sigma_prior[["df"]]
      sigma_scale <- solve(sigma_prior[["scale"]])
      # Drawing from the prior needs the prior to be proper, which a Wishart on
      # k variables only is from k degrees of freedom on. Fewer is a perfectly
      # good prior to estimate with -- the posterior adds one degree of freedom
      # per observation -- so this is a restriction on 'method', not on the
      # prior, and rWishart's own message does not say which of the two is at
      # fault.
      if (sigma_df < k) {
        stop("Initial values drawn from the prior need a proper prior for Sigma: ",
             "argument 'sigma$df' of add_priors() is ", sigma_df,
             " and must be at least the number of endogenous variables, ", k,
             ". Use method = \"maxlik\" or raise 'sigma$df'.")
      }
      object[["initial"]][["u_sigma_inv"]] <- matrix(stats::rWishart(1, df = sigma_df, Sigma = sigma_scale)[,,1], k)
    }
    
    if (error == "ald") {
      object <- .add_initial_values_ald(object, u = u, from_prior = TRUE)
    }
  }
  
  return(object)
}

# Initial values of the two blocks a quantile regression model has and a mean
# regression does not: the scale of the asymmetric Laplace, one per equation,
# and the latent scales of the mixture it is written as, one per observation
# and equation.
.add_initial_values_ald <- function(object, u, from_prior) {
  
  k <- object[["model"]][["k"]]
  tt <- ncol(u)
  
  if (from_prior) {
    scale_prior <- object[["priors"]][["u_scale"]]
    u_scale <- 1 / stats::rgamma(k, shape = scale_prior[["shape"]], rate = scale_prior[["rate"]])
  } else {
    # Given residuals, the maximum likelihood scale of an asymmetric Laplace is
    # the mean of the check function, which is what makes it the starting value
    # to take from a fit. Residuals of exactly zero would put it at zero, which
    # the sampler divides by, so it is floored.
    q <- object[["model"]][["quantile"]]
    u_scale <- apply(u, 1, function(z) {mean(z * (q - as.numeric(z < 0)))})
  }
  
  object[["initial"]][["u_scale"]] <- matrix(pmax(u_scale, 1e-8), k)
  # The latent scales are redrawn in the first sweep, before anything reads them
  # for their own sake, so they only have to be positive: they are the variance
  # the first draw of the coefficients is weighted by.
  object[["initial"]][["w"]] <- matrix(1, tt, k)
  
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
