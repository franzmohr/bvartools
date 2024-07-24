

# Extracts the posterior draws of the measurement error draws
.extract_posterior_sigma <- function(object) {
  
  draws_sigma <- NULL
  
  draws_sigma[["coeffs"]] <- t(object[["posteriors"]][["sigma"]][["coeffs"]])
  
  # if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
  #   k <- object[["model"]][["k"]]
  #   draws_c[["sigma"]] <- matrix(NA, n_c, draws)
  #   draws_c[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_c[1:n_c]])
  # }
  
  if ("lambda" %in% names(object[["posteriors"]][["sigma"]])) {
    k <- object[["model"]][["k"]]
    sigma_lambda <- matrix(diag(NA_real_, k), k * k, object[["model"]][["iterations"]])
    sigma_lambda[which(lower.tri(diag(1, k))), ] <- object[["posteriors"]][["sigma"]][["lambda"]]
    sigma_lambda[which(upper.tri(diag(1, k))), ] <- object[["posteriors"]][["sigma"]][["lambda"]]
    draws_sigma[["lambda"]] = sigma_lambda
  }
  
  return(draws_sigma)
}