#' Impulse Response Function
#' 
#' Computes the impulse response coefficients of an object of class "bvar" for
#' \code{n.ahead} steps.
#' 
#' @param object an object of class "bvar", usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param shock size of the impulse shock.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param type type of the impulse resoponse. Possible choices are forecast error "feir"
#' (default), orthogonalised "oir", structural "sir", and generalised "gir" impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' 
#' @details The function produces different types of impulse responses for the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' Forecast error impulse responses \eqn{\Phi_i} are obtained by recursions
#' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i-j} A_j,   i = 1, 2,...,}
#' with \eqn{\Phi_0 = I_K}.
#' 
#' Orthogonalised impulse responses \eqn{\Theta^o_i} are calculated as \eqn{\Theta^o_i = \Phi_i P,}
#' where P is the lower triangular Choleski decomposition of \eqn{\Sigma}.
#' 
#' Structural impulse responses \eqn{\Theta^s_i} are calculated as \eqn{\Theta^s_i = \Phi_i A_0^{-1}}.
#' 
#' Generalised impulse responses for variable \eqn{j}, i.e. \eqn{\Theta^g_ji} are calculated as
#' \eqn{\Theta^g_ji = \sigma_{jj}^{-1/2} \Phi_i \Sigma e_j}, where \eqn{\sigma_{jj}} is the variance
#' of the \eqn{j^{th}} diagonal element of \eqn{\Sigma} and \eqn{e_i} is a selection vector containing
#' one in its \eqn{j^{th}} element and zero otherwise.
#' 
#' @return A time-series object of class "bvarirf".
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
irf <- function(object, impulse = NULL, response = NULL, n.ahead = 5,
                shock = NULL, ci = .95, type = "feir", cumulative = FALSE) {
  
  if (!"bvar" %in% class(object)) {
    stop("Object must be of class 'bvar'.")
  }
  
  if (is.null(object$y) | is.null(dimnames(object$y)[[1]])) {
    stop("The argument 'object' must include a named matrix of endogenous variables.")
  }
  
  if (type %in% c("oir", "gir")) {
    if (is.null(object$Sigma)) {
      stop("OIR or GIR require to specify the 'Sigma' argument.")
    }
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
  }
  
  if (type == "sir") {
    if (is.null(object$A0)) {
      stop("Structural impulse response requires to specify the 'A0' argument.")
    }
    need_A0 <- TRUE
  } else {
    need_A0 <- FALSE
  }
  
  impulse <- which(dimnames(object$y)[[1]] == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(dimnames(object$y)[[1]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- NROW(object$y)
  store <- nrow(object$A)
  
  A <- c()
  for (i in 1:store) {
    temp <- list(A = matrix(object$A[i, ], k))
    if (need_Sigma) {
      temp$Sigma <- matrix(object$Sigma[i, ], k) 
    }
    if (need_A0) {
      temp$A0 = matrix(object$A0[i, ], k) 
    }
    A[[i]] <- temp
  }

  result <- lapply(A, .ir, h = n.ahead, type = type, impulse = impulse, response = response)
  
  result <- t(matrix(unlist(result), n.ahead + 1))
  
  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  pr <- c(ci_low, .5, ci_high)
  result <- stats::ts(t(apply(result, 2, stats::quantile, probs = pr)))
  stats::tsp(result) <- c(0, n.ahead, stats::tsp(result)[3])
  
  class(result) <- append("bvarirf", class(result))
  return(result)
}