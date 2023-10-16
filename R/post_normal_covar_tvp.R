#' Posterior Simulation of Error Covariance Coefficients
#' 
#' Produces posterior draws of time varying error covariance coefficients.
#' 
#' @param y a \eqn{K \times T} matrix of data with \eqn{K} as the number of
#' endogenous variables and \eqn{T} the number of observations.
#' @param u_omega_i matrix of error variances of the measurement equation.
#' Either a \eqn{K \times K} matrix for constant variances or
#' a \eqn{KT \times KT} matrix for time varying variances.
#' @param v_sigma_i matrix of error variances of the state equation.
#' Either an \eqn{M \times M} matrix for constant variances or
#' an \eqn{MT \times MT} matrix for time varying variances, where \eqn{M} is the
#' number of estimated variables.
#' @param psi_init a vector of inital values of the state equation.
#' 
#' @details For the multivariate model \eqn{A_{0,t} y_t = u_t} with \eqn{u_t \sim N(0, \Omega_t)}
#' the function produces a draw of the lower triangular part of \eqn{A_{0,t}} similar as in
#' Primiceri (2005), i.e., using \deqn{y_t = Z_t \psi_t + u_t,}
#' where 
#' \deqn{Z_{t} = \begin{bmatrix} 0 & \dotsm & \dotsm & 0 \\ -y_{1, t} & 0 & \dotsm & 0 \\ 0 & -y_{[1,2], t} & \ddots & \vdots \\ \vdots & \ddots & \ddots & 0 \\ 0 & \dotsm & 0 & -y_{[1,...,K-1], t} \end{bmatrix}}
#' and \eqn{y_{[1,...,K-1], t}} denotes the first to \eqn{(K-1)}th elements of the vector \eqn{y_t}.
#' 
#' The algorithm of Chan and Jeliazkov (2009) is used to obtain time varying coefficients.
#' 
#' @references
#' 
#' Chan, J., & Jeliazkov, I. (2009). Efficient simulation and integrated likelihood estimation in state space
#' models. \emph{International Journal of Mathematical Modelling and Numerical Optimisation, 1}(1/2), 101–120.
#' \doi{10.1504/IJMMNO.2009.030090}
#' 
#' Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy.
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
#' v_sigma_i <- diag(1000, 3)
#' psi_init <- matrix(0, 3)
#' 
#' # Obtain posterior draw
#' post_normal_covar_tvp(y, u_omega_i, v_sigma_i, psi_init)
#'
#' @export
post_normal_covar_tvp <- function(y, u_omega_i, v_sigma_i, psi_init) {
  
  k <- nrow(y)
  if (k == 1L) {
    stop("Argument 'y' must contain at least two variables.")
  }
  n_covar <- k * (k - 1) / 2
  tt <- ncol(y)
  y <- matrix(y)
  
  # Generate z for lower triangular design
  z <- .prep_covar_data(y, k, tt, TRUE)
  
  # Get positions of values, on which variables in z are regressed
  pos_used <- rep(k * 0:(tt - 1), each = k - 1) + 2:k
  
  # Trim endogenous variables
  y <- matrix(y[pos_used, ])
  
  # Trim 
  if (NCOL(u_omega_i) == k & NCOL(u_omega_i) == k) {
    u_omega_i <- kronecker(Matrix::Diagonal(tt, 1), u_omega_i)
  }
  # Trim error variance matrix
  u_omega_i <- u_omega_i[pos_used, pos_used]
  
  # Draw coefficients
  hh <- Matrix::Diagonal(n_covar * tt, 1)
  diag(hh[-(1:n_covar), -(n_covar * (tt - 1) + 1:n_covar)]) <- -1
  
  if (NCOL(v_sigma_i) == n_covar & NCOL(v_sigma_i) == n_covar) {
    v_sigma_i <- kronecker(Matrix::Diagonal(tt, 1), v_sigma_i)
  }
  
  hh_temp <- t(hh) %*% v_sigma_i %*% hh
  x_temp <- t(z) %*% u_omega_i
  v_post <- hh_temp + x_temp %*% z
  mu_post <- solve(v_post, hh_temp %*% kronecker(matrix(1, tt), psi_init) + x_temp %*% y)
  
  return(matrix(mu_post + solve(chol(v_post), matrix(stats::rnorm(n_covar * tt)))))
}