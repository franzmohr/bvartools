

.add_priors_sv_helper <- function(object, sigma, k) {
  
  object[["priors"]][["sigma"]][["type"]] <- "sv"
  object[["priors"]][["sigma"]][["mu"]] <- matrix(sigma[["mu"]], k)
  object[["priors"]][["sigma"]][["v_i"]] <- diag(sigma[["v_i"]], k)
  object[["priors"]][["sigma"]][["shape"]] <- matrix(sigma[["shape"]], k)
  object[["priors"]][["sigma"]][["rate"]] <- matrix(sigma[["rate"]], k)
  object[["priors"]][["sigma"]][["state_variance"]] <- matrix(sigma[["state_variance"]], k)
  object[["priors"]][["sigma"]][["offset"]] <- matrix(sigma[["offset"]], k)
  
  return(object)
}
