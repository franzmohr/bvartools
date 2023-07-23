#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Exponentially Weighted Moving Average
//' 
//' Produces a draw of log-volatilities.
//' 
//' @param u a \eqn{T \times K} matrix containing the time series.
//' @param lambda a scalar, vector or matrix of decay factors.
//' 
//' @details Some text.
//' 
//' @return A matrix with the same dimensions as \code{u}.
//' 
// [[Rcpp::export]]
arma::mat ewma(arma::mat u, const double lambda) {
   
   const int tt = u.n_rows;
   const int k = u.n_cols;
   
   arma::vec lambdavec = arma::vec(tt, arma::fill::value(lambda));
   lambdavec = arma::cumprod(lambdavec);
   
   arma::mat result = u * 0;
   arma::mat u_temp;

   // Square differences
   u = arma::pow(u, 2);

   for (int i = 0; i < tt; i++) {
     u_temp = arma::reverse(u.submat(0, 0, i, k - 1));
     u_temp = u_temp.each_col() % lambdavec.subvec(0, i);
     result.row(i) = arma::sum(u_temp, 0);
   }

   return (1 - lambda) * result;
 }

/*** R
library(AER); data("MarkPound"); y <- matrix(MarkPound); plot.ts(y)
#data("us_macrodata"); y <- us_macrodata; plot(y)
plot.ts(ewma(y, lambda = 0.94))
***/