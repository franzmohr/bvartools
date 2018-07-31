#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Draw Cointegration Matrix with Minnesota Prior
//' 
//' Produces a draw of coefficients of the cointegration matrix with a Minnesota Prior.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param ect Matrix containing the time series of the explanatory variables.
//' @param alpha Loading matrix of Pi.
//' @param Sigma_i Inverse covariance matrix.
//' 
//' @return List of the matrix \eqn{\Pi} and the cointegration matrix \eqn{\beta}.
//' 
//' @references
//' Koop, G., Leon-Gonzalez, R., & Strachan, R. W. (2010). Efficient posterior simulation for cointegrated models with priors on the cointegration space. \emph{Econometric Reviews}, 29(2), 224--242.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List posterior_koop2010(arma::mat y, arma::mat x, arma::mat ect, int r, arma::mat beta, arma::mat Sigma_i, arma::mat G_i, double v_i, arma::mat P_i, arma::mat B_mu_prior, arma::mat B_V_i_prior) {
  int n = y.n_rows;
  int nect = ect.n_rows;
  int nx = x.n_rows;
  int nalpha = n * r;
  int nbeta = nect * r;
  int ntot = n * (r + nx);
  int t = y.n_cols;
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  arma::mat ect_x;
  if (r > 0) {
    ect_x = arma::zeros<arma::mat>(r + nx, t);
    ect_x.rows(0, r - 1) = arma::trans(beta) * ect;
    ect_x.rows(r, r + nx - 1) = x;
  } else {
    ect_x = x;
  }

  arma::mat mu_prior = arma::zeros<arma::mat>(ntot, 1);
  mu_prior.rows(nalpha, ntot - 1) = B_mu_prior;
  arma::mat V_prior = arma::zeros<arma::mat>(ntot, ntot);
  if (r > 0) {
    V_prior.submat(0, 0, nalpha - 1, nalpha - 1) = arma::kron(v_i * arma::trans(beta) * P_i * beta, G_i);
  }
  V_prior.submat(nalpha, nalpha, ntot - 1, ntot - 1) = B_V_i_prior;
  
  arma::mat alpha_V_post = arma::inv(arma::kron(ect_x * arma::trans(ect_x), Sigma_i) + V_prior);
  arma::vec alpha_mu_post = alpha_V_post * (arma::vectorise(Sigma_i * y * arma::trans(ect_x)) + mu_prior);
  arma::svd(U, s, V, alpha_V_post);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  arma::vec z = arma::randn<arma::vec>(ntot);
  arma::mat B = reshape(alpha_mu_post + A * z, n, r + nx);
  
  arma::mat alpha;
  arma::mat Alpha;
  arma::mat Beta;
  arma::mat Pi;
  arma::vec Beta_mu_post;
  arma::mat Beta_V_post;
  if (r > 0) {
    alpha = B.cols(0, r - 1);
    arma::svd(U, s, V, alpha);
    //arma::mat Alpha = arma::zeros<arma::mat>(n, r);
    Alpha = U.cols(0, r - 1) * arma::trans(V);
    
    Beta_V_post = arma::inv(arma::kron(arma::trans(Alpha) * G_i * Alpha, v_i * P_i) + arma::kron(arma::trans(Alpha) * Sigma_i * Alpha, ect * arma::trans(ect)));
    Beta_mu_post = Beta_V_post * vectorise(ect * arma::trans(y) * Sigma_i * Alpha);
    
    arma::svd(U, s, V, Beta_V_post);
    A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
    z = arma::randn<arma::vec>(nbeta);
    Beta = arma::reshape(Beta_mu_post + A * z, nect, r);
    
    arma::svd(U, s, V, Beta);
    alpha = Alpha * V * arma::diagmat(s) * arma::trans(V);
    beta = U.cols(0, r - 1) * arma::trans(V);
    Pi = alpha * arma::trans(beta);
  } else {
    alpha = arma::zeros<arma::mat>(1,1);
    beta= arma::zeros<arma::mat>(1,1);
    Pi = arma::zeros<arma::mat>(n, nect);
  }
  
  arma::mat B_final = B.cols(r, r + nx - 1);
  arma::mat y_dot = y - Pi * ect - B_final * x;

  Rcpp::List result = Rcpp::List::create(B_final, Pi, alpha, beta, y_dot);
  return result;
}