#' Impulse Response Functions
#' 
#' Produces a draw of an impulse respone function.
#' 
#' @param object an object of class "bvars", usually, a result of a call to `bvars`.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param shock size of the impulse shock.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the credible intervals. Defaults to 0.95.
#' @param type type of the impulse resoponse. Possible choices are forecast error `feir`
#' (default), orthogonalised `oir`, structural `sir`, and generalised `gir` impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' 
#' @return A time-series object of class "irf.bvars".
#' 
#' @export
irf <- function(object, impulse = NULL, response = NULL, n.ahead = 5,
                shock = "sd", ci = .95, type = "feir", cumulative = FALSE) {
  
  if (!"bvars" %in% class(object)) {
    stop("Object must be of class 'bvars'.")
  }
  
  if (is.null(object$y) | is.null(dimnames(object$y)[[1]])) {
    stop("The argument 'object' must include a named matrix.")
  }
  n <- NROW(object$y)
  store <- nrow(object$A)
  
  A <- c()
  for (i in 1:store) {
    A <- c(A, list(matrix(object$A[i, ], n)))
  }
  
  temp <- lapply(A, feir, h = n.ahead)
  rm(A)
  
  impulse <- which(dimnames(object$y)[[1]] == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(dimnames(object$y)[[1]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  if (type == "feir") {
    result <- lapply(temp, function(x, h, impulse, response, n) {x[response + n * 0:h, impulse]},
                     h = n.ahead, impulse = impulse, response = response, n = n)
  }
  
  if (type == "sir") {
    if (is.null(object$A0)) {
      stop("Structural impulse response requires to specify the 'A0' argument.")
    }
    
    for (i in 1:store) {
      temp[[i]] <- list(Phi = temp[[i]],
                        A0 = matrix(object$A0[i, ], n))
    }
  }
  
  if (type %in% c("oir", "gir")) {
    if (is.null(object$Sigma)) {
      stop("OIR or GIR require to specify the 'Sigma' argument.")
    }
    
    for (i in 1:store) {
      temp[[i]] <- list(Phi = temp[[i]],
                        Sigma = matrix(object$Sigma[i, ], n))
    }
  }
  
  if (type %in% c("oir", "sir")) {
    result <- lapply(temp, sir, impulse = impulse, response = response, type = type)
  }
  
  if (type == "gir") {
    result <- lapply(temp, gir, impulse = impulse, response = response, shock = shock)
  }
  
  result <- t(matrix(unlist(result), n.ahead + 1))
  
  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }
  
  ci <- 1 - ci
  ci_low <- ci / 2
  ci_high <- 1 - ci / 2
  result <- stats::ts(t(apply(result, 2, stats::quantile, probs = c(ci_low, .5, ci_high))))
  
  stats::tsp(result) <- c(0, n.ahead, stats::tsp(result)[3])
  
  class(result) <- append("irf.bvars", class(result))
  return(result)
}