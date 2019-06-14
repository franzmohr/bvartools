#' Forecast Error Variance Decomposition
#' 
#' Produces the forecast error variance decomposition of a Bayesian VAR model.
#' 
#' @param object an object of class \code{"bvar"}, usually, a result of a call to \code{\link{bvar}}
#' or \code{\link{bvec_to_bvar}}.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param type type of the impulse responses used to calculate forecast error variable decompositions.
#' Possible choices are orthogonalised \code{oir} (default) and generalised \code{gir} impulse responses.
#' @param normalise_gir logical. Should the GIR-based FEVD be normalised?
#' 
#' @details The function produces forecast error variance decompositions (FEVD) for the VAR model
#' \deqn{y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + A_0^{-1} u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' If the FEVD is based on the orthogonalised impulse resonse (OIR), the FEVD will be calculated as
#' \deqn{\omega^{OIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i P e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\Phi_i} is the forecast error impulse response for the \eqn{i}th period,
#' \eqn{P} is the lower triangular Choleski decomposition of the variance-covariance
#' matrix \eqn{\Sigma}, \eqn{e_j} is a selection vector for the response variable and
#' \eqn{e_k} a selection vector for the impulse variable.
#'
#' If \code{type = "sir"}, the structural FEVD will be
#' calculated as \deqn{\omega^{SIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma^{\frac{1}{2}} e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1}  \Sigma A_0^{-1\prime} \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#'
#' If \code{type = "gir"}, the generalised FEVD will be
#' calculated as \deqn{\omega^{GIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#' 
#' If \code{type = "sgir"}, the structural generalised FEVD will be
#' calculated as \deqn{\omega^{SGIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma A_0^{-1\prime} \Phi_i^{\prime} e_j )}}.
#' 
#' Since GIR-based FEVDs do not add up to unity, they can be normalised by setting \code{normalise_gir = TRUE}.
#' 
#' @return A time-series object of class \code{"bvarfevd"}.

#' @examples
#' data("e1")
#' e1 <- diff(log(e1))
#' 
#' data <- gen_var(e1, p = 2, deterministic = "const")
#' 
#' y <- data$Y[, 1:73]
#' x <- data$Z[, 1:73]
#' 
#' set.seed(1234567)
#' 
#' iter <- 500 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations should be much higher, e.g. 30000.
#' 
#' burnin <- 100 # Number of burn-in draws
#' store <- iter - burnin
#' 
#' t <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' m <- k * nrow(x) # Number of estimated coefficients
#' 
#' # Set (uninformative) priors
#' a_mu_prior <- matrix(0, m) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, m) # Inverse of the prior covariance matrix
#' 
#' u_sigma_df_prior <- 0 # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
#' u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom
#' 
#' # Initial values
#' u_sigma_i <- diag(.00001, k)
#' u_sigma <- solve(u_sigma_i)
#' 
#' # Data containers for posterior draws
#' draws_a <- matrix(NA, m, store)
#' draws_sigma <- matrix(NA, k^2, store)
#' 
#' # Start Gibbs sampler
#' for (draw in 1:iter) {
#'   # Draw conditional mean parameters
#'   a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
#' 
#' # Draw variance-covariance matrix
#' u <- y - matrix(a, k) %*% x # Obtain residuals
#' u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
#' u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#' u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
#' 
#' # Store draws
#' if (draw > burnin) {
#'   draws_a[, draw - burnin] <- a
#'   draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#' 
#' # Generate bvar object
#' bvar_est <- bvar(y = y, x = x, A = draws_a[1:18,],
#'                  C = draws_a[19:21, ], Sigma = draws_sigma)
#' 
#' # Generate forecast error variance decomposition (FEVD)
#' bvar_fevd <- fevd(bvar_est, response = "cons")
#' 
#' # Plot FEVD
#' plot(bvar_fevd, main = "FEVD of consumption")

#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
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
  if (!type %in% c("oir", "sir", "gir", "sgir")) {
    stop("The specified type of the used impulse response is not known.")
  }
  if(is.null(response)) {
    stop("Please provide a valid response variable.")
  }
  
  k <- nrow(object$y)

  if (type %in% c("sgir", "sir")) {
    if (is.null(object$A0)) {
      stop("Structural impulse responses require that draws of 'A0' are contained in the 'bvar' object.")
    }
    A0_i <- solve(matrix(colMeans(object$A0), k))
  } else {
    A0_i <- diag(1, k)
  }
    
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
  
  if (type == "oir") {
    P <- t(chol(Sigma))
  }
  if (type == "sir") {
    P <- A0_i %*% sqrt(Sigma)
  }
  if (type %in% c("gir", "sgir")) {
    P <- A0_i %*% Sigma
  }
  
  numerator <- matrix(NA, k, n.ahead + 1)
  numerator[, 1] <- (ej_t %*% phi[1:k, ] %*% P)^2
  for (i in 2:(n.ahead + 1)) {
    numerator[, i] <- numerator[, i - 1] + (ej_t %*% phi[(i - 1) * k + 1:k, ] %*% P)^2
  }
  
  if (type == "gir") {
    numerator <- numerator / Sigma[response, response]
  }
  
  P <- NULL
  if (type %in% c("gir", "oir")) {
    P <- Sigma
  }
  if (type %in% c("sgir", "sir")) {
    P <- A0_i %*% tcrossprod(Sigma, A0_i)
  }
  
  mse <- rep(NA, n.ahead + 1)
  mse[1] <- ej_t %*% phi[1:k,] %*% P %*% t(phi[1:k,]) %*% t(ej_t)
  for (i in 2:(n.ahead + 1)) {
    mse[i] <- mse[i - 1] + ej_t %*% phi[(i - 1) * k + 1:k,] %*% P %*% t(phi[(i - 1) * k + 1:k,]) %*% t(ej_t)
  }

  result <- apply(numerator, 1, function(x, y) {x / y}, y = mse)
  
  if (type %in% c("gir", "sgir")) {
    if (normalise_gir) {
      result <- t(apply(result, 1, function(x) {x / sum(x)}))
    }
  }
  
  result <- stats::ts(result, start = 0, frequency = 1)
  dimnames(result) <- list(NULL, dimnames(object$y)[[1]])
  
  class(result) <- append("bvarfevd", class(result))
  return(result)
}