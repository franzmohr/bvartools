#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

//' Algorithm of Chan and Jeliazkov (2009)
//' 
//' A simplified implementation of the algorithm of Chan and Jeliazkov (2009).
//' 
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param z a sparse \eqn{KT \times MT} matrix of explanatory variables.
//' @param sigma_i a sparse \eqn{KT \times KT} inverse variance-covariance matrix of the measuresment
//' equation. For constant matrices \eqn{\Sigma^{-1}} the matrix corresponds to
//' \eqn{I_{T} \times \Sigma^{-1}}. Otherwise, the matrix is block diagonal.
//' @param b a sparse \eqn{MT \times MT} block diagonal difference matrix of the state equation.
//' @param q_i a sparse diagonal \eqn{MT \times MT} inverse variance-covariance matrix of the state
//' equation. For constant matrices \eqn{Q^{-1}} the matrix corresponds to
//' \eqn{I_{T} \times Q^{-1}}.
//' @param a0 an M-dimensional vector of initial states.
//' 
//' @details The function uses a simplified version of the algorithm of of Chan and Jeliazkov (2009)
//' to produce a draw of the state vector \eqn{a_t} for \eqn{t = 1,...,T} for a state space model
//' with measurement equation
//' \deqn{y_t = Z_t a_t + u_t}
//' and transition equation 
//' \deqn{a_{t} = a_{t - 1} + v_t,}
//' where \eqn{u_t \sim N(0, \Sigma_{t})} and \eqn{v_t \sim N(0, Q_{t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = z_t^{\prime} \otimes I_K} is a \eqn{K \times M} matrix of regressors with
//' \eqn{z_t} as a vector of regressors.
//' 
//' The implementation follows the depiction in chapter 20 of Chan, Koop, Poirier and Tobias (2019).
//' 
//' @return A \eqn{M \times T} matrix of state vector draws.
//' 
//' @examples
//' 
//' library(Matrix) # For sparse matrices
//' 
//' data("us_macrodata")
//' us_macrodata <- us_macrodata[time(us_macrodata) < 1990, ]
//' us_macrodata <- ts(us_macrodata, start = c(1959, 2), frequency = 4)
//' 
//' temp <- gen_var(us_macrodata, p = 2)
//' y <- temp$Y
//' x <- temp$Z
//' 
//' k <- NROW(y)
//' tt <- NCOL(y)
//' m <- NROW(x) * k
//' 
//' z <- matrix(0, k * tt, m * tt)
//' for (i in 1:tt) {
//'   z[(i - 1) * k + 1:k, (i - 1) * m + 1:m] <- kronecker(t(x[, i]), diag(1, k))
//' }
//' z <- Matrix(z, sparse = TRUE)
//' 
//' # Initial state
//' b0 <- tcrossprod(y, x) %*% solve(tcrossprod(x))
//' u <- y - b0 %*% x
//' b0 <- matrix(b0)
//' 
//' # Initial variance-covariance matrix
//' sigma_i <- solve(tcrossprod(u) / tt)
//' sigma_i <- kronecker(diag(1, tt), sigma_i)
//' sigma_i <- Matrix(sigma_i, sparse = TRUE)
//' 
//' # Initial variances of the state equation
//' q_i <- diag(1 / .001, m)
//' q_i <- kronecker(diag(1, tt), q_i)
//' q_i <- Matrix(q_i, sparse = TRUE)
//' 
//' # Differences matrix of state equation
//' h <- diag(1, tt * m)
//' diag(h[-(1:m), -(((tt - 1) * m):(tt * m))]) <- -1
//' h <- Matrix(h, sparse = TRUE)
//' 
//' # Draw states
//' b <- chan_jeliazkov(y = y, z = z, sigma_i = sigma_i,
//'                     b = h, q_i = q_i, a0 = b0)
//' 
//' @references
//' 
//' Chan, J. C. C., & Jeliazkov, I. (2009). Efficient simulation and integrated likelihood
//' estimation in state space models. \emph{International Journal of Mathematical Modelling
//' and Numerical Optimisation, 1}(1/2), 101--120. \url{https://doi.org/10.1504/ijmmno.2009.030090}
//' 
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
// [[Rcpp::export]]
Rcpp::List chan_jeliazkov(arma::mat y, arma::sp_mat z, arma::sp_mat sigma_i,
                         arma::sp_mat b, arma::sp_mat q_i, arma::vec a0) {
  
  arma::uword k = y.n_rows;
  int tt = y.n_cols;
  int nvars = z.n_cols / tt;
  
  y = arma::reshape(y, k * tt, 1);
  arma::sp_mat HQiH = arma::trans(b) * q_i * b;
  arma::sp_mat Zsi = arma::trans(z) * sigma_i;
  arma::sp_mat K = HQiH + Zsi * z;
  arma::vec mu = arma::spsolve(K, HQiH * arma::kron(arma::ones<arma::vec>(tt), a0) + Zsi * y, "lapack");
  arma::mat Kd(K);
  arma::sp_mat A = arma::sp_mat(arma::trans(arma::chol(Kd, "lower")));
  arma::mat a = mu + arma::spsolve(A, arma::randn<arma::vec>(nvars * tt), "lapack");
  arma::mat u = arma::reshape(y - z * a, k, tt);
  
  return Rcpp::List::create(Rcpp::Named("a") = arma::join_rows(a0, arma::reshape(a, nvars, tt)),
                            Rcpp::Named("u") = u);
}