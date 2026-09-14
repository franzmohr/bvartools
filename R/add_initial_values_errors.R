
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
    # The samplers read shape and rate as a Gamma(shape, rate) prior on each
    # error precision, so that is what is drawn. This used to draw the inverse
    # of a gamma with half the shape and the inverse rate, which for shape = 50
    # and rate = 25 started the precisions near 0.001 rather than near 2.
    if (error %in% c("gamma", "gamma+covar")) {
      sigma_shape <- as.numeric(sigma_prior[["shape"]])
      if (any(sigma_shape <= 0)) {
        stop("Initial values drawn from the prior need a proper gamma prior on the error ",
             "precisions: argument 'sigma$shape' of add_priors() must be larger than 0.",
             call. = FALSE)
      }
      object[["initial"]][["u_omega_inv"]] <- diag(stats::rgamma(k, shape = sigma_shape,
                                                                 rate = as.numeric(sigma_prior[["rate"]])), k)
    }

    if (error %in% c("sv", "sv+covar")) {
      mu <- sigma_prior[["mu"]]
      vinv <- sigma_prior[["v_inv"]]
      h_draw <- .draw_normal_prior(mu, vinv, "sigma")
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

# One draw from a normal prior given by its mean and precision matrix, which is
# how add_priors() stores it.
#
# With R'R = V^-1, the draw mu + R^-1 e has covariance (R'R)^-1 = V. It used to
# be mu + R e, whose covariance is the precision: a prior standard deviation of
# 0.1 started the chain at coefficients with a standard deviation of 10. A prior
# that is uninformative for some element has nothing to draw from, and chol()
# only said that a leading minor was not positive.
.draw_normal_prior <- function(mu, v_inv, argument) {
  if (any(diag(as.matrix(v_inv)) <= 0)) {
    stop("Initial values cannot be drawn from a prior that is uninformative for some ",
         "of its elements: the prior precision in argument '", argument, "' of ",
         "add_priors() must be larger than 0 for every element. Use the least squares ",
         "initial values instead.", call. = FALSE)
  }
  mu + backsolve(chol(v_inv), stats::rnorm(length(mu)))
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
#
# The state precisions start at the mean of their gamma priors unless the
# initial values are to be drawn from the prior. They used to be drawn under
# every method, and that draw was the only use of the RNG in the "ols" and
# "maxlik" paths: a seed set between add_initial_values() and
# add_posterior_coefficients() -- the natural place for it -- then fixed the
# sampler but not the point its chain started from, and a TVP model did not
# repeat its draws across R sessions.
.add_initial_values_state_errors <- function(object, method) {

  # The samplers read the prior of each precision as Gamma(shape, rate), whose
  # mean is shape / rate. The draw used to be from Gamma(shape / 2, rate / 2),
  # which has the same mean and twice the variance.
  state_precision <- function(prior) {
    shape <- as.numeric(prior[["shape"]])
    rate <- as.numeric(prior[["rate"]])
    n <- length(shape)
    result <- diag(0, n)
    for (i in 1:n) {
      if (method == "prior") {
        result[i, i] <- stats::rgamma(1, shape = shape[i], rate = rate[i])
      } else {
        result[i, i] <- shape[i] / rate[i]
      }
    }
    result
  }

  if (object[["model"]][["tvp"]] & !is.null(object[["data"]][["train"]][["z"]])) {
    object[["initial"]][["a_sigma_inv"]] <- state_precision(object[["priors"]][["a"]])
  }

  use_covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  if (object[["model"]][["tvp"]] & use_covar & object[["model"]][["k"]] > 1) {
    object[["initial"]][["psi_sigma_inv"]] <- state_precision(object[["priors"]][["psi"]])
  }

  return(object)
}
