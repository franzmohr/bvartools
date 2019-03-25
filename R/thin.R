#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bvar"} or \code{"bvec"}.
#' 
#' @param object an object of class \code{"bvar"} or \code{"bvec"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' 
#' @return An object of class \code{"bvar"} or \code{"bvec"}.
#' 
#' @export
thin <- function(object, thin = 5) {
  if (!any(class(object) %in% c("bvar", "bvec"))) {
    stop("Argument 'object' must be of class 'bvar' or 'bvec'.")
  }
  
  if (any(class(object) == "bvar")) {
    var <- TRUE
  }
  if (any(class(object) == "bvec")) {
    var <- FALSE
  }
  
  if (var) {
    draws <- nrow(object$A)
  } else {
    draws <- nrow(object$Pi)
  }
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  if (var) {
    object$A <- coda::mcmc(object$A[pos_thin,], start = start, end = end, thin = thin)
    if(!is.null(object$A0)) {
      object$A0 <- coda::mcmc(object$A0[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$B)) {
      object$B <- coda::mcmc(object$B[pos_thin,], start = start, end = end, thin = thin)
    }
  } else {
    object$Pi <- coda::mcmc(object$Pi[pos_thin,], start = start, end = end, thin = thin)
    if(!is.null(object$A0)) {
      object$A0 <- coda::mcmc(object$A0[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$alpha)) {
      object$alpha <- coda::mcmc(object$alpha[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$beta)) {
      object$beta <- coda::mcmc(object$beta[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Pi_x)) {
      object$Pi_x <- coda::mcmc(object$Pi_x[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Pi_d)) {
      object$Pi_d <- coda::mcmc(object$Pi_d[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Gamma)) {
      object$Gamma <- coda::mcmc(object$Gamma[pos_thin,], start = start, end = end, thin = thin)
    }
    if(!is.null(object$Upsilon)) {
      object$Upsilon <- coda::mcmc(object$Upsilon[pos_thin,], start = start, end = end, thin = thin)
    }
  }
  if(!is.null(object$C)) {
    object$C <- coda::mcmc(object$C[pos_thin,], start = start, end = end, thin = thin)
  }
  if(!is.null(object$Sigma)) {
    object$Sigma <- coda::mcmc(object$Sigma[pos_thin,], start = start, end = end, thin = thin)
  }
  return(object)
}