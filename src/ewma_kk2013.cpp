#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Exponentially Weighted Moving Average
//' 
//' Produces a matrix of covariance matrices based on the EWMA algorithm used
//' in Koop and Korobilis (2013).
//' 
//' @param u a \eqn{T \times k} matrix containing the time series.
//' @param lambda a numeric specifying the decay factor.
//' @param sigmainit a matrix with the initial values of sigma.
//' 
//' @details The 
//' \deqn{\Sigma_{t} = \lambda \Sigma_{t - 1} + (1 - \lambda) u_t u_t^{\prime}.}
//' 
//' @return A matrix.
//' 
//' @references
//' 
//' Koop, G., & Korobilis, D. (2013). Large time-varying parameter VARs.
//' \emph{Journal of Econometrics 177}(2), 185--198. \doi{10.1016/j.jeconom.2013.04.007}
//' 
// [[Rcpp::export]]
arma::mat ewma_kk2013(arma::mat u, double lambda, arma::mat sigmainit) {
   
   const int tt = u.n_rows;
   const int k = u.n_cols;
   arma::mat sigma = arma::zeros<arma::mat>(tt * k, k);
   
   sigma.submat(0, 0, k - 1, k - 1) = lambda * sigmainit + (1 - lambda) * u.row(0).t() * u.row(0); 
   for (int i = 1; i < tt; i++) {
     sigma.submat(i * k, 0, (i + 1) * k - 1, k - 1) = lambda * sigma.submat((i - 1) * k, 0, i * k - 1, k - 1) + (1 - lambda) * u.row(i).t() * u.row(i);
   }
   
   return sigma;
 }

/*** R
library(AER); data("MarkPound"); y <- matrix(MarkPound)
#data("us_macrodata"); y <- us_macrodata
plot.ts(y)
res <- ewma_kk2013(u = y,
                   lambda = 0.94,
                   sigmainit = var(y))
plot.ts(res)
# 
# head(res)
# sigma1 <- 0.94 * var(us_macrodata) + (1 - 0.94) * tcrossprod(us_macrodata[1, ])
# 0.94 * sigma1 + (1 - 0.94) * tcrossprod(us_macrodata[2, ])
***/