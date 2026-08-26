#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.log_likelihood_normal)]]
arma::mat log_likelihood_normal(const int k, const arma::mat u, const arma::mat sigma_inv) {
  
  
  int draws = sigma_inv.n_rows;
  int u_nrows = u.n_rows;
  int tt = u_nrows / k;
  int sigma_ncols = sigma_inv.n_cols;
  bool tvp_errors = sigma_ncols > (k * k);
  
  double part_a = -k * log(2 * arma::datum::pi) / 2 ;
  double part_b, part_c;
  const arma::mat diag_k = arma::eye(k, k);
  arma::mat u_sigma_inv;
    
  arma::mat loglik = arma::zeros<arma::mat>(draws, tt);
  
  for (int draw = 0; draw < draws; draw++)
  {
    if (tvp_errors) {
      for (int i = 0; i < tt; i++)
      {
        u_sigma_inv = arma::reshape(sigma_inv.submat(draw, i * k, draw, (i + 1) * k - 1), k, k);
        part_b = -log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
        part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
        loglik(draw, i) = part_a + part_b + part_c;
      }
    } else {
      u_sigma_inv = arma::reshape(sigma_inv.row(draw), k, k);
      part_b = -log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
      for (int i = 0; i < tt; i++)
      {
        part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
        loglik(draw, i) = part_a + part_b + part_c;
      }
    }
  }
  
  return loglik;
}