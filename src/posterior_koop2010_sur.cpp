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
// [[Rcpp::export]]
Rcpp::List posterior_koop2010_sur(arma::mat y, arma::mat x, arma::mat ect, arma::mat beta,
                                  arma::mat Sigma_i, arma::mat G_i, double v_i, arma::mat P_i,
                                  arma::mat B_mu_prior, arma::mat B_V_i_prior) {
  int n = y.n_rows;
  int t = y.n_cols;
  int nx = x.n_rows;
  int nect = ect.n_rows;
  int r = beta.n_cols;
  int nalpha = n * r;
  int nbeta = nect * r;
  int ntot = n * (nx + r);
  bool const_var = true;
  arma::mat S_i = arma::zeros<arma::mat>(n * t, n * t);
  if (Sigma_i.n_rows > n) {
    const_var = false;
  }
  
  arma::mat ect_x = arma::zeros<arma::mat>(r + nx, t);
  ect_x.rows(0, r - 1) = arma::trans(beta) * ect;
  ect_x.rows(r, r + nx - 1) = x;
  
  arma::mat Z = arma::kron(arma::trans(ect_x), arma::eye<arma::mat>(n, n));
  arma::mat Alpha = arma::zeros<arma::mat>(n, r);
  
  arma::mat ZHi = arma::zeros<arma::mat>(ntot, n * t);
  arma::mat ZHZ = arma::zeros<arma::mat>(ntot, ntot);
  arma::mat ZHy = arma::zeros<arma::mat>(ntot, 1);
  
  arma::mat U;
  arma::vec s;
  arma::mat V;
  
  if (const_var) {
    ZHi = arma::trans(Z) * arma::kron(arma::eye<arma::mat>(t, t), Sigma_i);
  } else {
    for (int i = 0; i < t; i++){
      S_i.submat(i * n, i * n, (i + 1) * n - 1, (i + 1) * n - 1) = Sigma_i.rows(i * n, (i + 1) * n - 1);
    }
    ZHi = arma::trans(Z) * S_i;
  }
  
  ZHZ = ZHi * Z;
  ZHy = ZHi * arma::vectorise(y);
  
  arma::mat mu_prior = arma::zeros<arma::mat>(ntot, 1);
  arma::mat V_prior = arma::zeros<arma::mat>(ntot, ntot);
  mu_prior.rows(nalpha, ntot - 1) = B_mu_prior;
  V_prior.submat(0, 0, nalpha - 1, nalpha - 1) = arma::kron(v_i * arma::trans(beta) * P_i * beta, G_i);
  V_prior.submat(nalpha, nalpha, ntot - 1, ntot - 1) = B_V_i_prior;

  arma::mat V_post = arma::inv(ZHZ + V_prior);
  arma::vec mu_post = V_post * (ZHy + mu_prior);
  
  arma::svd(U, s, V, V_post);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  arma::vec z = arma::randn<arma::vec>(ntot);
  arma::mat B = reshape(mu_post + A * z, n, r + nx);
  
  arma::mat alpha = B.cols(0, r - 1);
  arma::mat B_final = B.cols(r, r + nx - 1);
  
  arma::svd(U, s, V, alpha);
  Alpha = U.cols(0, r - 1) * arma::trans(V);
  
  ZHZ = arma::zeros<arma::mat>(nbeta, nbeta);
  ZHy = arma::zeros<arma::mat>(nbeta, 1);
  if (const_var) {
    for (int i = 0; i < t; i++){
      Z = arma::kron(Alpha, arma::trans(ect.col(i)));
      ZHZ = ZHZ + arma::trans(Z) * Sigma_i * Z;
      ZHy = ZHy + arma::trans(Z) * Sigma_i * y.col(i);
    } 
  } else {
    for (int i = 0; i < t; i++){
      Z = arma::kron(Alpha, arma::trans(ect.col(i)));
      ZHZ = ZHZ + arma::trans(Z) * Sigma_i.rows(i * n, (i + 1) * n - 1) * Z;
      ZHy = ZHy + arma::trans(Z) * Sigma_i.rows(i * n, (i + 1) * n - 1) * y.col(i);
    }
  }

  arma::mat Beta_V_post = arma::inv(ZHZ + arma::kron(arma::trans(Alpha) * G_i * Alpha, v_i * P_i));
  arma::vec Beta_mu_post = Beta_V_post * ZHy;

  arma::svd(U, s, V, Beta_V_post);
  A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  z = arma::randn<arma::vec>(nbeta);
  arma::mat Beta = arma::reshape(Beta_mu_post + A * z, nect, r);
  
  arma::svd(U, s, V, Beta);
  alpha = Alpha * V * arma::diagmat(s) * arma::trans(V);
  beta = U.cols(0, r - 1) * arma::trans(V);
  
  arma::mat Pi = alpha * arma::trans(beta);
  Rcpp::List result = Rcpp::List::create(B_final, Pi, alpha, beta);
  return result;
}