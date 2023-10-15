// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export(.prep_covar_data)]]
arma::sp_mat prep_covar_data(arma::vec y, int k, int tt, bool tvp) {
  
  int n_covar = k * (k - 1) / 2;
  arma::sp_mat z;
  
  if (tvp) {
    z = arma::zeros<arma::sp_mat>((k - 1) * tt, n_covar * tt);
    for (int j = 0; j < tt; j++) {
      for (int i = 1; i < k; i++) {
        z.submat(j * (k - 1) + i - 1,
                 j * n_covar + i * (i - 1) / 2,
                 j * (k - 1) + i - 1,
                 j * n_covar + (i + 1) * i / 2 - 1) = -arma::trans(y.subvec(j * k, j * k + i - 1));
      }
    }
  } else {
    z = arma::zeros<arma::sp_mat>((k - 1) * tt, n_covar);
    for (int j = 0; j < tt; j++) {
      for (int i = 1; i < k; i++) {
        z.submat(j * (k - 1) + i - 1,
                 i * (i - 1) / 2,
                 j * (k - 1) + i - 1,
                 (i + 1) * i / 2 - 1) = -arma::trans(y.subvec(j * k, j * k + i - 1));
      }
    }
  }
  
  return z;
}

/*** R
.prep_covar_data(matrix(1:9), 3, 3, TRUE)
*/