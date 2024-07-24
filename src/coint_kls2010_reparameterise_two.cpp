#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Cointegration Reparameterisation
//' 
//' Performs the second transformation of the loading and
//' cointegration matrix as proposed in Koop et al. (2010).
//' 
//' @param alpha a \eqn{K \times r} matrix.
//' @param beta an \eqn{M \times r} matrix.
//' 
//' @details The function performs two transformations:
//' \itemize{
//'   \item \eqn{A = \alpha (\beta^{\prime} \beta)^{1/2}}
//'   \item \eqn{B = \beta (\beta^{\prime} \beta)^{-1/2}}
//' }
//' 
//' @references
//' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
//' simulation for cointegrated models with priors on the cointegration space.
//' \emph{Econometric Reviews, 29}(2), 224-242. \doi{10.1080/07474930903382208}
//' 
//' @return A list of two matrices:
//' \item{alpha}{Loading matrix \eqn{A}.}
//' \item{beta}{Cointegration matrix \eqn{B} (semiorthogonal).}
//' 
//' @examples
//' 
//' # Generate input data
//' alpha <- matrix(c(-0.07, 0.17), 2)
//' beta <- matrix(c(1, -4), 2)
//' 
//' # Reparameterise
//' coint_kls2010_reparameterise_two(alpha, beta)
//' 
// [[Rcpp::export(coint_kls2010_reparameterise_two)]]
Rcpp::List coint_kls2010_reparameterise_two(const arma::mat alpha, const arma::mat beta) {
  
  const arma::mat BB_sqrt = arma::sqrtmat_sympd(arma::trans(beta) * beta);
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha * BB_sqrt,
                            Rcpp::Named("beta") = beta * arma::solve(BB_sqrt, arma::eye<arma::mat>(alpha.n_cols, alpha.n_cols)));
}

/*** R

alpha <- matrix(c(-0.07, 0.17), 2)
beta <- matrix(c(1, -4), 2)

coint_kls2010_reparameterise_two(alpha, beta)

*/