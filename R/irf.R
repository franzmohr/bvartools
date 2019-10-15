#' Impulse Response Function
#' 
#' Computes the impulse response coefficients of an object of class \code{"bvar"} for
#' \code{n.ahead} steps.
#' 
#' @param object an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param type type of the impulse resoponse. Possible choices are forecast error \code{"feir"}
#' (default), orthogonalised \code{"oir"}, structural \code{"sir"}, generalised \code{"gir"},
#' and structural generalised \code{"sgir"} impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the posterior impulse response function. Defaults to \code{FALSE} so that
#' the median and the credible intervals of the posterior draws are returned.
#' 
#' @details The function produces different types of impulse responses for the VAR model
#' \deqn{y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + A_0^{-1} u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' Forecast error impulse responses \eqn{\Phi_i} are obtained by recursions
#' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i-j} A_j,   i = 1, 2,...,h}
#' with \eqn{\Phi_0 = I_K}.
#' 
#' Orthogonalised impulse responses \eqn{\Theta^o_i} are calculated as \eqn{\Theta^o_i = \Phi_i P},
#' where P is the lower triangular Choleski decomposition of \eqn{\Sigma}. \eqn{A_0} is assumed to
#' be an identity matrix.
#' 
#' Structural impulse responses \eqn{\Theta^s_i} are calculated as \eqn{\Theta^s_i = \Phi_i A_0^{-1}}.
#' 
#' (Structural) Generalised impulse responses for variable \eqn{j}, i.e. \eqn{\Theta^g_ji} are calculated as
#' \eqn{\Theta^g_{ji} = \sigma_{jj}^{-1/2} \Phi_i A_0^{-1} \Sigma e_j}, where \eqn{\sigma_{jj}} is the variance
#' of the \eqn{j^{th}} diagonal element of \eqn{\Sigma} and \eqn{e_i} is a selection vector containing
#' one in its \eqn{j^{th}} element and zero otherwise. If the \code{"bvar"} object does not contain draws
#' of \eqn{A_0}, it is assumed to be an identity matrix.
#' 
#' @return A time-series object of class \code{"bvarirf"} and if \code{keep_draws = TRUE} a simple matrix.
#' 
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
#' # Generate impulse response
#' IR <- irf(bvar_est, impulse = "income", response = "cons", n.ahead = 8)
#' 
#' # Plot
#' plot(IR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
#'
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
irf <- function(object, impulse = NULL, response = NULL, n.ahead = 5,
                ci = .95, type = "feir", cumulative = FALSE, keep_draws = FALSE) {
  
  if (!type %in% c("feir", "oir", "gir", "sir", "sgir")) {
    stop("Argument 'type' not known.")
  }
  
  if (!"bvar" %in% class(object)) {
    stop("Object must be of class 'bvar'.")
  }
  
  if (is.null(object$y) | is.null(dimnames(object$y)[[1]])) {
    stop("Argument 'object' must include a named matrix of endogenous variables.")
  }

  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (is.null(object$A0)) {
      stop("Structural impulse responses require that draws of 'A0' are contained in the 'bvar' object.")
    }
    need_A0 <- TRUE
  }
  
  if (type %in% c("oir", "gir", "sgir")) {
    if (is.null(object$Sigma)) {
      stop("OIR, GIR, SGIR require that draws of 'Sigma' are contained in the 'bvar' object.")
    }
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
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
  
  if (!is.null(attr(object$y, "ts_info"))) {
    freq <- attr(object$y, "ts_info")[3]
  } else {
    freq <- 1
  }
  
  if (!keep_draws) {
    ci_low <- (1 - ci) / 2
    ci_high <- 1 - ci_low
    pr <- c(ci_low, .5, ci_high)
    result <- stats::ts(t(apply(result, 2, stats::quantile, probs = pr)), start = 0, frequency = freq) 
  }
  
  class(result) <- append("bvarirf", class(result))
  return(result)
}