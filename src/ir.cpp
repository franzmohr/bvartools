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
  
  arma::mat Sigma, temp;
  arma::mat P = arma::eye<arma::mat>(k, k);
  arma::mat norm_error;
  
  if (type == fe) {
    P = P * Rcpp::as<double>(A["shock"]);
  }
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"])));
    P = P * arma::diagmat(1 / arma::vectorise(P.diag())) * Rcpp::as<double>(A["shock"]);
  }
  if (type == sir) {
    P = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k)) * Rcpp::as<double>(A["shock"]);
  }
  if (type == gir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = Sigma / arma::as_scalar(Sigma(impulse - 1, impulse - 1)) * Rcpp::as<double>(A["shock"]);
  }
  if (type == sgir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k));
    P = P * Sigma / arma::as_scalar(Sigma(impulse - 1, impulse - 1)) * Rcpp::as<double>(A["shock"]);
  }

  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp;
  arma::vec theta = arma::zeros<arma::vec>(h + 1);
  
  // Initial value of theta
  temp = phi_temp * P;
  theta(0) = arma::as_scalar(temp(response - 1, impulse - 1));
  
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = coef.cols(0, p - 1);
  
  for (int i = 1; i <= h; i++) {
    // FEIR
    phi_temp.zeros();
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp;
    
    // Potential transformation of FEIR
    temp = phi_temp * P;
    theta(i) = arma::as_scalar(temp(response - 1, impulse - 1));
  }

  return theta;
}