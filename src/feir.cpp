#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Forecast Error Impulse Response
//' 
//' Produces a forecast error impulse response matrix of a vector autoregressive model.
//' 
//' @param A \eqn{n x np} matrix of coefficients, where \eqn{n} is the number of endogenous
//' variables and \eqn{p} is the number of lags.
//' @param h integer specifying the steps.
//' 
//' @details The function produces the moving average representation \eqn{\Phi_i} for the model
//' \deqn{y_t = \sum_{i}^{p} A_i y_{t - i} + \epsilon_{t}}
//' by recursions
//' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i - j} A_j}
//' for \eqn{i = 1, 2, ...} and starting with \eqn{\Phi_0 = I_K}.
//' 
//' @return a \eqn{nh x n} matrix.
//' 
// [[Rcpp::export]]
arma::mat feir(arma::mat A, int h = 5) {
  int n = A.n_rows;
  int p = A.n_cols / n;
  
  if (h < p) {
    p = h * n;
  } else {
    p = A.n_cols;
  }
  
  arma::mat A_temp = arma::zeros<arma::mat>(n, h * n);
  A_temp.cols(0, p - 1) = A.cols(0, p - 1); 
  
  // Generate output matrix
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * n, n);
  phi.rows(0, n - 1) = arma::eye<arma::mat>(n, n);
  
  arma::mat phi_temp = arma::zeros<arma::mat>(n, n);
  
  for (int i = 1; i <= h; i++) {
    phi_temp.zeros();
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * n, (i - j + 1) * n - 1) * A_temp.cols((j - 1) * n, j * n - 1);
    }
    phi.rows(i * n, (i + 1) * n - 1) = phi_temp;
  }
  
  return phi;
}