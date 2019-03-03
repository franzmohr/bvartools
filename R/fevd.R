#' Forecast Error Variance Decomposition
#' 
#' Produces the forecast error variance decomposition of a Bayesian VAR model.
#' 
#' @param object an object of class "bvar", usually, a result of a call to \code{\link{bvar}}
#' or \code{\link{bvec_to_bvar}}.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param type type of the impulse responses used to calculate forecast error variable decompositions.
#' Possible choices are orthogonalised \code{oir} (default) and generalised \code{gir} impulse responses.
#' @param normalise_gir logical. Should the GIR-based FEVD be normalised?
#' 
#' @details The function produces forecast error variance decompositions (FEVD) for the VAR model
#' \deqn{y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' If the FEVD is based on the orthogonalised impulse resonse (OIR), the FEVD will be calculated as
#' \deqn{\omega^{OIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i P e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\Phi_i} is the forecast error impulse response for the \eqn{i}th period,
#' \eqn{P} is the lower triangular Choleski decomposition of the variance-covariance
#' matrix \eqn{\Sigma}, \eqn{e_j} is a selection vector for the response variable and
#' \eqn{e_k} a selection vector for the impulse variable.
#'
#' If \code{type = "gir"}, the FEVD will be
#' calculated as \deqn{\omega^{GIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#' 
#' Since GIR-based FEVDs do not add up to unity, they can be normalised by setting \code{normalise_gir = TRUE}.
#' 
#' @return A time-series object of class "bvarfevd".
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
fevd <- function(object, response = NULL, n.ahead = 5, type = "oir", normalise_gir = FALSE) {
  if (!"bvar" %in% class(object)) {
    stop("Object must be of class 'bvar'.")
  }
  if (is.null(object$y) | is.null(dimnames(object$y)[[1]])) {
    stop("The argument 'object' must include a named matrix of endogenous variables.")
  }
  if (is.null(object$Sigma)) {
    stop("The 'bvar' object must include draws of the variance-covariance matrix Sigma.")
  }
  if (!type %in% c("oir", "gir")) {
    stop("The specified type of the used impulse response is not known.")
  }
  if(is.null(response)) {
    stop("Please provide a valid response variable.")
  }
  
  k <- nrow(object$y)
  
  A_temp <- matrix(colMeans(object$A), k)
  A <- matrix(0, k, k * n.ahead)
  A[, 1:ncol(A_temp)] <- A_temp
  phi <- matrix(0, k * (1 + ncol(A) / k), k)
  phi[1:k, 1:k] <- diag(1, k)
  for (i in 1:(n.ahead)) {
    phi_temp <- matrix(0, k, k)
    for (j in 1:i) {
      phi_temp = phi_temp + phi[(i - j) * k + 1:k,] %*% A[, (j - 1) * k + 1:k];
    }
    phi[i * k + 1:k,] <- phi_temp
  }
  Sigma <- matrix(colMeans(object$Sigma), k)
  
  response <- which(dimnames(object$y)[[1]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  ej_t <- matrix(0, 1, k)
  ej_t[,response] <- 1
  
  mse <- rep(NA, n.ahead + 1)
  mse[1] <- ej_t %*% phi[1:k,] %*% Sigma %*% t(phi[1:k,]) %*% t(ej_t)
  for (i in 2:(n.ahead + 1)) {
    mse[i] <- mse[i - 1] + ej_t %*% phi[(i - 1) * k + 1:k,] %*% Sigma %*% t(phi[(i - 1) * k + 1:k,]) %*% t(ej_t)
  }
  
  if (type == "oir") {
    P <- t(chol(Sigma))
  }
  if (type == "gir") {
    P <- Sigma
  }
  
  numerator <- matrix(NA, k, n.ahead + 1)
  numerator[, 1] <- (ej_t %*% phi[1:k, ] %*% P)^2
  for (i in 2:(n.ahead + 1)) {
    numerator[, i] <- numerator[, i - 1] + (ej_t %*% phi[(i - 1) * k + 1:k, ] %*% P)^2
  }

  result <- apply(numerator, 1, function(x, y) {x / y}, y = mse)
  if (type == "gir") {
    result <- result / Sigma[response, response]
    if (normalise_gir) {
      result <- t(apply(result, 1, function(x) {x / sum(x)}))
    }
  }
  dimnames(result) <- list(NULL, dimnames(object$y)[[2]])
  result <- stats::ts(result, start = 0, frequency = 1)
  
  class(result) <- append("bvarfevd", class(result))
  return(result)
}