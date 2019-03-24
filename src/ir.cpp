#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.ir)]]
arma::vec ir(Rcpp::List A, int h, std::string type, int impulse, int response) {
  arma::mat coef = Rcpp::as<arma::mat>(A["A"]);
  
  int k = coef.n_rows;
  int p = coef.n_cols / k;
  
  std::string fe ("feir");
  std::string oir ("oir");
  std::string sir ("sir");
  std::string gir ("gir");
  std::string sgir ("sgir");
  
  if (h < p) {
    p = h * k;
  } else {
    p = coef.n_cols;
  }
  
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = coef.cols(0, p - 1); 
  
  arma::mat P, Sigma, temp;
  
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp;
  arma::vec theta = arma::zeros<arma::vec>(h + 1);
  
  if (type == fe) {
    theta(0) = arma::as_scalar(phi(response - 1, impulse - 1));
  }
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"])));
  }
  if (type == sir) {
    P = arma::inv(Rcpp::as<arma::mat>(A["A0"]));
  }
  if (type == gir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = Sigma / sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1)));
  }
  if (type == sgir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = arma::inv(Rcpp::as<arma::mat>(A["A0"])) * Sigma / sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1)));
  }
  if (type != fe) {
    temp = phi_temp * P;
    theta(0) = arma::as_scalar(temp(response - 1, impulse - 1));
  }
  
  for (int i = 1; i <= h; i++) {
    phi_temp.zeros();
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp;
    
    if (type == fe) {
      theta(i) = arma::as_scalar(phi_temp(response - 1, impulse - 1));
    } else {
      temp = phi_temp * P;
      theta(i) = arma::as_scalar(temp(response - 1, impulse - 1));
    }
  }

  return theta;
}