#include <RcppArmadillo.h>

#include "impact_matrix.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.ir)]]
arma::vec ir(Rcpp::List A, int h, std::string type, int impulse, int response) {
  
  arma::mat coef = Rcpp::as<arma::mat>(A["A"]);
  
  int k = coef.n_rows;
  const int lags = coef.n_cols / k;
  
  std::string fe ("feir");
  std::string oir ("oir");
  std::string sir ("sir");
  std::string gir ("gir");
  std::string sgir ("sgir");
  std::string custom ("custom");

  if (type != fe && type != oir && type != sir && type != gir &&
      type != sgir && type != custom) {
    Rcpp::stop("ir: unknown type \"%s\".", type);
  }

  // Regressor columns the recursion can reach: lags beyond the horizon never
  // enter it. This is a count of columns, not a lag order, and it is zero on
  // impact -- where there is no recursion and nothing to slice.
  int n_use = coef.n_cols;
  if (h < lags) {
    n_use = h * k;
  }
  
  arma::mat Sigma, temp;
  arma::mat P = arma::eye<arma::mat>(k, k);
  arma::mat norm_error;
  
  if (type == fe) {
    P = P * Rcpp::as<double>(A["shock"]);
  }
  // The orthogonalised and generalised types count a shock in standard
  // deviations: a column of the Choleski factor, and Sigma e_j over the
  // standard deviation of the impulse variable, are the responses to a shock of
  // one standard deviation. That is Phi_i P and sigma_jj^{-1/2} Phi_i Sigma e_j
  // as the documentation states them, and the convention of vars::irf.
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"]))) * Rcpp::as<double>(A["shock"]);
  }
  if (type == sir) {
    P = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k)) * Rcpp::as<double>(A["shock"]);
  }
  if (type == gir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = Sigma / std::sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1))) * Rcpp::as<double>(A["shock"]);
  }
  if (type == sgir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k));
    P = P * Sigma / std::sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1))) * Rcpp::as<double>(A["shock"]);
  }
  if (type == custom) {
    // The caller has already decided what a shock of size one is, so `shock`
    // only rescales it. Nothing here normalises the columns of P the way the
    // "oir" branch above does -- an identification that means to deliver unit
    // shocks has to arrive that way.
    P = bvartools_impact_matrix(A, k, "irf") * Rcpp::as<double>(A["shock"]);
  }

  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp;
  arma::vec theta = arma::zeros<arma::vec>(h + 1);
  
  // Initial value of theta
  temp = phi_temp * P;
  theta(0) = arma::as_scalar(temp(response - 1, impulse - 1));
  
  // On impact the response is Phi_0 P = P, which theta(0) already holds, so the
  // coefficients are only needed once the recursion below actually runs.
  arma::mat A_temp;
  if (h > 0) {
    A_temp = arma::zeros<arma::mat>(k, h * k);
    if (n_use > 0) {
      A_temp.cols(0, n_use - 1) = coef.cols(0, n_use - 1);
    }
  }
  
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