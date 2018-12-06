#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

//' Durbin and Koopman Simulation Smoother
//' 
//' The function is an implementation of the Kalman filter and backward smoothing
//' algorithm proposed by Durbin and Koopman (2002).
//' 
//' @param y a \eqn{n x T} matrix containing the time series of the dependent variable.
//' @param Z a \eqn{nT x m} matrix containing the time series of the explanatory variables.
//' @param Sigma_v a \eqn{n x n} or \eqn{nT x n} variance-covariance matrix.
//' @param Sigma_w a \eqn{m x m} variance-covariance matrix of the transition equation.
//' @param B a \eqn{m x m} autocorrelation function of the transition equation.
//' @param a_init a \eqn{m x 1} vector of the initial values.
//' @param Sigma_w_init a \eqn{m x m} variance-covariance matrix of the initial parameter values.
//' 
//' @details For the state space model with the measurement equation
//' \deqn{y_t = Z_t a_t + v_t}
//' and the transition equation 
//' \deqn{a_{t+1} = B_t a_t + w_t,}
//' where \eqn{v_t \sim N(0, \Sigma_{v_t})} and \eqn{w_t \sim N(0, \Sigma_{w_t})} the
//' function produces a draw of the state vector \eqn{a_t} for \eqn{T = 1,...,T}.
//' 
//' @return A \eqn{m x T} matrix of parameter values.
//' 
//' @references
//' 
//' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for state space time series analysis. \emph{Biometrika}, 89(3), 603--615.
//' 
// [[Rcpp::export]]
arma::mat kalman_dk(arma::mat y, arma::mat Z, arma::mat Sigma_v, arma::mat Sigma_w, arma::mat B, arma::vec a_init, arma::mat Sigma_w_init) {
  // Algorithm 2  of Durbin and Koopman (2002)
  int n = y.n_rows;
  int t = y.n_cols;
  int nvars = Z.n_cols;
  
  arma::mat H_temp = arma::zeros<arma::mat>(n * t, n);
  if (Sigma_v.n_rows == n){
    for (int i = 0; i < t; i++){
      H_temp.rows(i * n, (i + 1) * n - 1) = Sigma_v;
    }
    Sigma_v = H_temp;
  }
  
  arma::mat Q_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (Sigma_w.n_rows == nvars){
    for (int i = 0; i < t; i++){
      Q_temp.rows(i * nvars, (i + 1) * nvars - 1) = Sigma_w;
    }
    Sigma_w = Q_temp;
  }
  
  arma::mat T_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (B.n_rows == nvars){
    for (int i = 0; i < t; i++){
      T_temp.rows(i * nvars, (i + 1) * nvars - 1) = B;
    }
    B = T_temp;
  }
  
  arma::mat yplus = arma::zeros<arma::mat>(n, t);
  arma::mat aplus = arma::zeros<arma::mat>(nvars, t + 1);
  
  arma::mat U;
  arma::mat V;
  arma::vec s;
  
  svd(U, s, V, Sigma_w_init);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  arma::vec z = arma::randn<arma::vec>(nvars);
  aplus.col(0) = a_init + A * z;
  
  int p1 = 0;
  int p2 = 0;
  int pA1 = 0;
  int pA2 = 0;
  for (int i = 0; i < t; i++){
    p1 = i * n;
    p2 = (i + 1) * n - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    
    svd(U, s, V, Sigma_w.rows(pA1, pA2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
    z = arma::randn<arma::vec>(nvars);
    aplus.col(i + 1) = B.rows(pA1, pA2) * aplus.col(i) + A * z;
    
    svd(U, s, V, Sigma_v.rows(p1, p2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
    z = arma::randn<arma::vec>(n);
    yplus.col(i) =  Z.rows(p1, p2) * aplus.col(i) + A * z;
  }
  arma::mat ystar = y - yplus;
  
  arma::mat a = arma::zeros<arma::mat>(nvars, t + 1);
  a.col(0) = a_init;
  arma::mat P = Sigma_w_init;
  arma::mat v = arma::zeros<arma::mat>(n, t);
  arma::mat Fi = arma::zeros<arma::mat>(n * t, n);
  arma::mat K = arma::zeros<arma::mat>(nvars * t, n);
  arma::mat L = arma::zeros<arma::mat>(nvars * t, nvars);
  for (int i = 0; i < t ; i++){
    p1 = i * n;
    p2 = (i + 1) * n - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    v.col(i) = ystar.col(i) - Z.rows(p1, p2) * a.col(i);
    Fi.rows(p1, p2) = inv(Z.rows(p1, p2) * P * arma::trans(Z.rows(p1, p2)) + Sigma_v.rows(p1, p2));
    K.rows(pA1, pA2) = B.rows(pA1, pA2) * P * arma::trans(Z.rows(p1, p2)) * Fi.rows(p1, p2);
    L.rows(pA1, pA2) = B.rows(pA1, pA2) - K.rows(pA1, pA2) * Z.rows(p1, p2);
    a.col(i + 1) = B.rows(pA1, pA2) * a.col(i) + K.rows(pA1, pA2) * v.col(i);
    P = B.rows(pA1, pA2) * P * arma::trans(L.rows(pA1, pA2)) + Sigma_w.rows(pA1, pA2);
  }
  
  arma::mat r = arma::zeros<arma::mat>(nvars, t);
  for (int i = (t - 1); i > 0; i--){
    r.col(i - 1) = arma::trans(Z.rows(i * n, (i + 1) * n - 1)) * Fi.rows(i * n, (i + 1) * n - 1) * v.col(i) + arma::trans(L.rows(i * nvars, (i + 1) * nvars - 1)) * r.col(i);
  }
  arma::vec r0 = arma::trans(Z.rows(0, n - 1)) * Fi.rows(0, n - 1) * v.col(0) + arma::trans(L.rows(0, nvars - 1)) * r.col(0);
  
  arma::mat ahatstar = arma::zeros<arma::mat>(nvars, t + 1);
  ahatstar.col(0) = a_init + Sigma_w_init * r0;
  for (int i = 0; i < t; i++){
    ahatstar.col(i + 1) = B.rows(i * nvars, (i + 1) * nvars - 1) * ahatstar.col(i) + Sigma_w.rows(i * nvars, (i + 1) * nvars - 1) * r.col(i);
  }
  return ahatstar + aplus;
}