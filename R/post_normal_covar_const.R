#' Posterior Simulation of Error Covariance Coefficients
#' 
#' Produces posterior draws of constant error covariance coefficients.
#' 
#' @param y a \eqn{K \times T} matrix of data with \eqn{K} as the number of
#' endogenous variables and \eqn{T} the number of observations.
#' @param u_omega_i matrix of error variances of the measurement equation.
#' Either a \eqn{K \times K} matrix for constant variances or
#' a \eqn{KT \times KT} matrix for time varying variances.
#' @param prior_mean vector of prior means. In case of TVP, this vector is used
#' as initial condition.
#' @param prior_covariance_i inverse prior covariance matrix. In case of TVP, this matrix
#' is used as initial condition.
#' 
#' @details For the multivariate model \eqn{A_0 y_t = u_t} with \eqn{u_t \sim N(0, \Omega_t)}
#' the function produces a draw of the lower triangular part of \eqn{A_0} similar as in
#' Primiceri (2005), i.e., using \deqn{y_t = Z_t \psi + u_t,}
#' where 
#' \deqn{Z_{t} = \begin{bmatrix} 0 & \dotsm & \dotsm & 0 \\ -y_{1, t} & 0 & \dotsm & 0 \\ 0 & -y_{[1,2], t} & \ddots & \vdots \\ \vdots & \ddots & \ddots & 0 \\ 0 & \dotsm & 0 & -y_{[1,...,K-1], t} \end{bmatrix}}
#' and \eqn{y_{[1,...,K-1], t}} denotes the first to \eqn{(K-1)}th elements of the vector \eqn{y_t}.
#' 
#' @references Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy.
#' \emph{The Review of Economic Studies, 72}(3), 821--852. \doi{10.1111/j.1467-937X.2005.00353.x}
#' 
#' @returns A matrix.
#' 
#' @examples
#' # Load example data
#' data("e1")
#' y <- log(t(e1))
#' 
#' # Generate artificial draws of other matrices
#' u_omega_i <- diag(1, 3)
#' prior_mean <- matrix(0, 3)
#' prior_covariance_i <- diag(0, 3)
#' 
#' # Obtain posterior draw
#' post_normal_covar_const(y, u_omega_i, prior_mean, prior_covariance_i)
#'
#' @export
post_normal_covar_const <- function(y, u_omega_i, prior_mean, prior_covariance_i) {
  
  k <- nrow(y)
  if (k == 1L) {
    stop("Argument 'y' must contain at least two variables.")
  }
  n_covar <- k * (k - 1) / 2
  tt <- ncol(y)
  y <- matrix(y)
  
  # Generate z for lower triangular design
  # Use C++ to speed up the for loop
  z <- .prep_covar_data(y, k, tt, FALSE)
  
  # Get positions of values, on which variables in z are regressed
  pos_used <- rep(k * 0:(tt - 1), each = k - 1) + 2:k
  
  # Trim endogenous variables
  y <- matrix(y[pos_used, ])
  
  # Adapt error variance matrix
  if (NCOL(u_omega_i) == k & NCOL(u_omega_i) == k) {
    u_omega_i <- kronecker(Diagonal(tt, 1), u_omega_i)
  }
  # Trim error variance matrix
  u_omega_i <- u_omega_i[pos_used, pos_used]
  
  # Draw coefficients
  v_post <- prior_covariance_i + crossprod(z, u_omega_i) %*% z
  mu_post <- solve(v_post, prior_covariance_i %*% prior_mean + crossprod(z, u_omega_i) %*% y)
  
  return(matrix(mu_post + solve(chol(v_post), matrix(rnorm(n_covar)))))
}