#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Posterior Data Preparation
//' 
//' Generates a lower triangular block matrix with ones on the main diagonal,
//' where the off-diagonal elements are blocks.
//' 
//' @param a \eqn{K^2p}-dimensional vector of coefficients. See 'Details'.
//' @param k integer of the number of columns per block.
//' @param tt integer specifying 
//' 
//' @return A sparse matrix.
//' 
//' @details For the \eqn{K \times Kp} matrix A, where \eqn{a = vec(A)}, with
//' \eqn{A = \left[A_1, A_2, ..., A_p\right]} the function constructs the
//' following sparse \eqn{KT \times KT} diagonal block matrix:
//' \deqn{\begin{bmatrix} I_{K}  & 0 & 0 & 0 & \dots & 0 \\-A_{1} & I_{K} & 0 & 0 & \dots & 0 \\-A_{2} & -A_{1}& I_{K} & 0 & \dots & 0 \\ 0 & -A_{2}& -A_{1}& I_{K} & \dots & 0 \\ \vdots & \ddots& \ddots& \ddots& \ddots& 0 \\  0 & \dots & 0 & -A_{2}& -A_{1}& I_{K} \end{bmatrix}.}
//' 
//' @examples
//' a <- matrix(1:8)
//' generate_lower_block_diagonal(a, 2, 5)
//' 
// [[Rcpp::export(generate_lower_block_diagonal)]]
arma::sp_mat generate_lower_block_diagonal(arma::mat& a, int& k, int& tt) {

  if (a.n_elem % (k * k) != 0) {
    Rcpp::stop("Length of vector 'a' is not a multiple of k * k.");
  }
  
  int p = a.n_elem / (k * k);
  // Transform matrix
  arma::mat amat = arma::reshape(a, k, k * p);
  arma::mat amat_t = arma::zeros<arma::mat>(k * p, k);
  for (int i = 0; i < p; i++) {
    amat_t.rows(i * k, (i + 1) * k - 1) = -amat.cols(i * k, (i + 1) * k - 1); // Minus sign!
  }

  // Update H_a
  arma::sp_mat result = arma::speye<arma::sp_mat>(tt * k, tt * k);
  for (int i = 0; i < (tt - 1); i++) {
    if (i < tt - p) {
      result.submat((i + 1) * k, i * k,
                    (i + 1) * k + k * p - 1, (i + 1) * k - 1) = amat_t;
    } else {
      result.submat((i + 1) * k, i * k,
                    (i + 1) * k + k * (p - (tt - i - 1)) - 1, (i + 1) * k - 1) = amat_t.rows(0, k * (p - (tt - i - 1)) - 1);
    }
  }

  return result;
}
