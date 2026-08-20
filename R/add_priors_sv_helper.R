

.add_priors_sv_helper <- function(object, sigma, k) {
  
  object[["priors"]][["u_sigma"]][["type"]] <- "sv"
  object[["priors"]][["u_sigma"]][["mu"]] <- matrix(sigma[["mu"]], k)
  object[["priors"]][["u_sigma"]][["v_inv"]] <- diag(sigma[["v_i"]], k)
  object[["priors"]][["u_sigma"]][["shape"]] <- matrix(sigma[["shape"]], k)
  object[["priors"]][["u_sigma"]][["rate"]] <- matrix(sigma[["rate"]], k)
  object[["priors"]][["u_sigma"]][["sigma"]] <- matrix(sigma[["state_variance"]], k)
  object[["priors"]][["u_sigma"]][["offset"]] <- matrix(sigma[["offset"]], k)
  
  return(object)
}
