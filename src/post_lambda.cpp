#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.post_lambda)]]
arma::mat post_lambda(arma::mat& x, arma::mat& ff, arma::mat& prior_vinv, arma::sp_mat& uinv, arma::mat lambda) {

  int m = x.n_cols;
  int n = ff.n_rows;

  int nlambda_i, pos1, pos2;
  arma::mat K_l, lambda_hat, lambda_temp;

  // Draw lambda
  int lambda_count = 0;
  for (int i = 1; i < m; i++) { // For each row in lambda, where coefficient should be estimated
    pos1 = lambda_count;
    if (i <= (n - 1)) { // If number of estimated coefficients is lower than n_factor
      nlambda_i = i;
      pos2 = lambda_count + i - 1;
      K_l = prior_vinv.submat(pos1, pos1, pos2, pos2) + ff.rows(0, i - 1) * ff.rows(0, i - 1).t() * uinv(i, i);
      lambda_hat = arma::solve(K_l, ff.rows(0, i - 1) * (x.col(i) - arma::trans(ff.row(i))) * uinv(i, i));
    } else {
      pos2 = lambda_count + n - 1;
      nlambda_i = n;
      K_l = prior_vinv.submat(pos1, pos1, pos2, pos2) + ff * ff.t() * uinv(i, i);
      lambda_hat = arma::solve(K_l, ff * x.col(i) * uinv(i, i));
    }
    lambda_temp = arma::trans(lambda_hat + arma::solve(chol(K_l, "lower"), arma::randn(nlambda_i)));
    lambda.submat(i, 0, i, nlambda_i - 1) = lambda_temp;
    lambda_count = pos2 + 1;
  }

  return lambda;

}
