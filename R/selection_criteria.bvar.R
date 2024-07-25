#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class \code{"bvar"}.
#'
#' @param object an object of class \code{"bvar"}, usually, the
#' output of a call to \code{\link{draw_posterior}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list.
#' 
#' @export
selection_criteria.bvar <- function(object, ...){
  
  result <- list()
  
  # Get model specs
  tvp <- any(unlist(object[["specifications"]][["tvp"]]))
  sv <- object[["specifications"]][["tvp"]][["Sigma"]]
  tt <- NROW(object[["y"]])
  k <- object[["specifications"]][["dims"]][["K"]]
  p <- object[["specifications"]][["lags"]][["p"]]
  result[["p"]] <- object[["specifications"]][["lags"]][["p"]]
  if (object[["specifications"]][["dims"]][["M"]] != 0) {
    result[["s"]] <- object[["specifications"]][["lags"]][["s"]] 
  }
  
  if (tvp) {
    temp_pars <- list()
    length(temp_pars) <- tt
  } else {
    temp_pars <- NULL
  }
  
  if (sv) {
    draws <- nrow(object[["Sigma"]][[1]])
  } else {
    draws <- nrow(object[["Sigma"]])
  }
  
  x <- NULL
  if (!is.null(object[["x"]])) {
    
    vars <- c("A", "B", "C")
    for (j in vars) {
      if (!is.null(object[[j]])) {
        if (is.list(object[[j]])) {
          for (period in 1:tt) {
            temp_pars[[period]] <- cbind(temp_pars[[period]], object[[j]][[period]]) 
          }
        } else {
          temp_pars <- cbind(temp_pars, object[[j]]) 
        }
      }
    } 
    
    z <- kronecker(object[["x"]], diag(1, k))
    if (tvp) {
      for (period in 1:tt) {
        temp_pars[[period]] <- t(temp_pars[[period]]) 
      }
      temp_pars <- do.call("rbind", temp_pars)
      z <- sur_const_to_tvp(z, k, tt)
    } else {
      temp_pars <- t(temp_pars)
    }
    
    u <- t(matrix(matrix(t(object[["y"]])), k * tt, draws) - z %*% temp_pars)
    tot_pars <- NCOL(object[["x"]])
    
  } else {
    u <- t(matrix(matrix(t(object[["y"]])), k * tt, draws))
    tot_pars <- 0
  }
  
  LL <- matrix(NA, draws, tt) # Get LogLik
  if (sv) {
    sigma <- matrix(NA_real_, k * tt, k)
  } else {
    sigma <- matrix(NA_real_, k, k)
  }
  
  for (j in 1:draws) {
    
    if (sv) {
      for (period in 1:tt) {
        sigma[(period - 1) * k + 1:k,] <- matrix(object[["Sigma"]][[period]][j,], k)
      }
    } else {
      sigma <- matrix(object[["Sigma"]][j,], k)
    }
    
    # LogLik
    LL[j, ] <- loglik_normal(matrix(u[j, ], k, tt), sigma)
  }
  
  
  ll <- sum(colMeans(LL))
  
  result[["LL"]] <- ll
  result[["AIC"]] <- 2 * tot_pars - 2 * ll
  result[["BIC"]] <- log(tt) * tot_pars - 2 * ll
  result[["HQ"]] <- 2 * log(log(tt)) * tot_pars - 2 * ll
  
  result <- as.data.frame(result)
  
  result <- list(ll = LL,
                 npara = tot_pars,
                 summary = result)
  
  class(result) <- append("selcrit", class(result))
  
  return(result)
}