

# Extracts the posterior draws of the structural parameters
.extract_posterior_a0 <- function(object) {
  
  draws_a0 <- NULL
  
  if (object[["model"]][["structural"]]) {
    
    k <- object[["model"]][["k"]]
    tt <- nrow(object[["data"]][["y"]])
    draws <- object[["model"]][["iterations"]]
    n_z <- ncol(object[["data"]][["z"]])
    pos <- which(lower.tri(diag(1, k)))
    pos_a0 <- (n_z - length(pos) + 1):n_z
    
    if (object[["model"]][["tvp"]]) {
      pos_a0_long <- rep(pos_a0, tt) + rep(0:(tt - 1), each = length(pos_a0)) * n_z
      draws_a0[["coeffs"]] <- matrix(diag(1, k), k * k * tt, draws)
      draws_a0[["coeffs"]][rep(0:(tt - 1) * k * k, each = length(pos)) + rep(pos, tt), ] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_a0_long])
    } else {
      draws_a0[["coeffs"]] <- matrix(diag(1, k), k * k, draws)
      draws_a0[["coeffs"]][pos, ] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_a0])
    }
    
    if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
      draws_a0[["sigma"]] <- matrix(0, k * k, draws)
      draws_a0[["sigma"]][pos, ] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_a0])
    }
    
    if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
      draws_a0[["lambda"]] <- matrix(diag(1, k), k * k, draws)
      draws_a0[["lambda"]][pos, ] <- object[["posteriors"]][["a"]][["lambda"]][, pos_a0]
      draws_a0[["lambda"]][-pos, ] <- NA_real_
    }
  }
  
  return(draws_a0)
}