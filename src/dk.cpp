#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

//' Durbin and Koopman Kalman Filter
//' 
//' Produces a draw from the Kalman filter proposed by Durbin and Koopman.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param Z Matrix containing the time series of the explanatory variables.
//' @param H Matirx.
//' @param Q sfg.
//' @param T The inverse of the current draw of the variance-covariance matrix.
//' @param a1 fg.
//' @param P1 fg.
//' 
// [[Rcpp::export]]
arma::mat dk(arma::mat y, arma::mat Z, arma::mat H, arma::mat Q, arma::mat T, arma::vec a1, arma::mat P1) {
  // Algorithm2  of Durbin and Koopman (2002)
  int n = y.n_rows;
  int t = y.n_cols;
  int nvars = Z.n_cols;
  
  if (H.n_rows==n){
    arma::mat H_temp = arma::zeros<arma::mat>(n*t,n);
    for (int i = 0; i < t; i++){
      H_temp.rows(i * n, (i + 1) * n - 1) = H;
    }
    H = H_temp;
  }
  
  if (Q.n_rows == nvars){
    arma::mat Q_temp = arma::zeros<arma::mat>(nvars * t, nvars);
    for (int i = 0; i < t; i++){
      Q_temp.rows(i * nvars, (i + 1) * nvars - 1) = Q;
    }
    Q = Q_temp;
  }
  
  if (T.n_rows == nvars){
    arma::mat T_temp = arma::zeros<arma::mat>(nvars * t, nvars);
    for (int i = 0; i < t; i++){
      T_temp.rows(i * nvars, (i + 1) * nvars - 1) = T;
    }
    T = T_temp;
  }
  
  // Algorithm 2
  
  arma::mat yplus = arma::zeros<arma::mat>(n, t);
  arma::mat aplus = arma::zeros<arma::mat>(nvars, t + 1);

  arma::mat U;
  arma::mat V;
  arma::vec s;
  
  svd(U, s, V, P1);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  arma::vec z = arma::randn<arma::vec>(nvars);
  aplus.col(0) = a1 + A * z;
  
  int p1 = 0;
  int p2 = 0;
  int pA1 = 0;
  int pA2 = 0;
  for (int i = 0; i < t; i++){
    p1 = i * n;
    p2 = (i + 1) * n - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    
    svd(U, s, V, Q.rows(pA1, pA2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
    z = arma::randn<arma::vec>(nvars);
    aplus.col(i + 1) = T.rows(pA1, pA2) * aplus.col(i) + A * z;
    
    svd(U, s, V, H.rows(p1, p2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
    z = arma::randn<arma::vec>(n);
    yplus.col(i) =  Z.rows(p1, p2) * aplus.col(i) + A * z;
  }
  arma::mat ystar = y - yplus;
  
  arma::mat a = arma::zeros<arma::mat>(nvars, t + 1);
  a.col(0) = a1;
  arma::mat P = P1;
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
    //F.rows(p1, p2) = Z.rows(p1, p2) * P * arma::trans(Z.rows(p1, p2)) + H.rows(p1, p2);
    //Fi.rows(p1, p2) = inv(F.rows(p1, p2));
    Fi.rows(p1, p2) = inv(Z.rows(p1, p2) * P * arma::trans(Z.rows(p1, p2)) + H.rows(p1, p2));
    K.rows(pA1, pA2) = T.rows(pA1, pA2) * P * arma::trans(Z.rows(p1, p2)) * Fi.rows(p1, p2);
    L.rows(pA1, pA2) = T.rows(pA1, pA2) - K.rows(pA1, pA2) * Z.rows(p1, p2);
    a.col(i + 1) = T.rows(pA1, pA2) * a.col(i) + K.rows(pA1, pA2) * v.col(i);
    P = T.rows(pA1, pA2) * P * arma::trans(L.rows(pA1, pA2)) + Q.rows(pA1, pA2);
  }
  
  arma::mat r = arma::zeros<arma::mat>(nvars, t);
  for (int i = (t - 1); i > 0; i--){
    r.col(i - 1) = arma::trans(Z.rows(i * n, (i + 1) * n - 1)) * Fi.rows(i * n, (i + 1) * n - 1) * v.col(i) + arma::trans(L.rows(i * nvars, (i + 1) * nvars - 1)) * r.col(i);
  }
  arma::vec r0 = arma::trans(Z.rows(0, n - 1)) * Fi.rows(0, n - 1) * v.col(0) + arma::trans(L.rows(0, nvars - 1)) * r.col(0);
  
  arma::mat ahatstar = arma::zeros<arma::mat>(nvars, t + 1);
  ahatstar.col(0) = a1 + P1 * r0;
  for (int i = 0; i < t; i++){
    ahatstar.col(i + 1) = T.rows(i * nvars, (i + 1) * nvars - 1) * ahatstar.col(i) + Q.rows(i * nvars, (i + 1) * nvars - 1) * r.col(i);
  }
  arma::mat atilde = ahatstar + aplus;
  return atilde;
}