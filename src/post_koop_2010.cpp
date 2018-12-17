#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw of Cointegration Parameters
//' 
//' Produces a draw from a posterior density Koop et al. (2010).
//' 
//' @param y a \eqn{k \times T} matrix of differenced dependent variables.
//' @param beta a \eqn{k_{ect} \times r} matrix of variables in the error correction term.
//' @param ect a \eqn{k_{ect} \times T} matrix of regressor variables in levels and deterministic terms in the error correction term.
//' @param x a \eqn{k_{x} \times T} matrix of differenced regressor variables and unrestriced deterministic terms.
//' @param Sigma_i a \eqn{k \times k} variance-covariance matrix.
//' @param Gamma_mu_prior a \eqn{kk_{x} \times 1} prior mean vector of non-cointegration coefficients.
//' @param Gamma_V_i_prior a \eqn{kk_{x} \times kk_{x}} inverted prior covariance vector of non-cointegration coefficients.
//' @param v_i a numeric between 0 and 1 specifying the shrinkage.
//' @param P_tau_i a \eqn{k_{ect} \times k_{ect}} inverted matrix specifying the central location of \eqn{sp(\beta)}.
//' @param G_i a \eqn{k \times k} matrix.
//' 
//' @details The function produces a draw of coefficients for the model
//' \deqn{\Delta y_t = \alpha \beta' ECT_{t} + \Gamma X_{t} + u_t,}
//' where \eqn{\Delta y_t} is a \eqn{k \times 1} vector of differenced dependent variables,
//' \eqn{\alpha} is a \eqn{k \times r} loading matrix,
//' \eqn{\beta} is a \eqn{k_{ect} \times r} cointegration matrix,
//' \eqn{ECT_t} is a \eqn{k_{ect} \times 1} vector of regressors in the error correction term,
//' \eqn{\Gamma} is a \eqn{k \times k_{x}} matrix of non-cointegration parameters,
//' \eqn{X_t} is a \eqn{k_{x} \times 1} vector of non-cointegration regressors,
//' and \eqn{u_t} is a \eqn{k \times 1} error term.
//' 
//' @return A list containing the following elements:
//' \item{alpha}{a draw of the \eqn{k \times r} loading matrix.}
//' \item{beta}{a draw of the \eqn{k_{ect} \times r} cointegration matrix.}
//' \item{Pi}{a draw of the \eqn{k \times k_{ect}} matrix \eqn{\Pi = \alpha \beta'}.}
//' \item{Gamma}{a draw of the \eqn{k \times k_{x}} coefficient matrix for non-cointegration parameters.}
//' 
//' @references
//' 
//' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
//' simulation for cointegrated models with priors on the cointegration space.
//' \emph{Econometric Reviews}, 29(2), 224-242. \url{https://doi.org/10.1080/07474930903382208}
//' 
// [[Rcpp::export]]
Rcpp::List post_koop_2010(arma::mat y, arma::mat beta, arma::mat ect,
                          arma::mat x,
                          arma::mat Sigma_i,
                          arma::vec Gamma_mu_prior,
                          arma::mat Gamma_V_i_prior,
                          double v_i, arma::mat P_tau_i, arma::mat G_i){
  
  int k = y.n_rows;
  int r = beta.n_cols;
  int k_a = k * r;
  int k_b = ect.n_rows * r;
  
  int k_x = 0;
  int k_g = 0;
  bool incl_x = false;
  if (x.n_cols == y.n_cols) {
    incl_x = true;
    k_x = x.n_rows;
    k_g = k * k_x;
  }
  
  int k_ag = k_a + k_g;
  arma::mat Z = arma::zeros<arma::mat>(r + k_x, y.n_cols);
  Z.rows(0, r - 1) = arma::trans(beta) * ect;
  if (incl_x) {
    Z.rows(r, r + k_x - 1) = x; 
  }
  arma::mat S_ag_post = arma::kron(Z * arma::trans(Z), Sigma_i);
  arma::mat S_ag_prior = S_ag_post * 0;
  S_ag_prior.submat(0, 0, k_a - 1, k_a - 1) = arma::kron(v_i * (arma::trans(beta) * P_tau_i * beta), G_i);
  if (incl_x) {
    S_ag_prior.submat(k_a, k_a, k_ag - 1, k_ag - 1) = Gamma_V_i_prior; 
  }
  S_ag_post = arma::inv(S_ag_post + S_ag_prior);
  arma::vec mu_ag_post = S_ag_post * arma::reshape(Sigma_i * y * arma::trans(Z), k_ag, 1);
  
  arma::mat U_ag;
  arma::mat V_ag;
  arma::vec s_ag;
  arma::svd(U_ag, s_ag, V_ag, S_ag_post);
  arma::mat ag_sqrt = U_ag * arma::diagmat(sqrt(s_ag)) * arma::trans(V_ag);
  arma::vec z_ag = arma::randn<arma::vec>(k_ag);
  arma::mat ag = mu_ag_post + ag_sqrt * z_ag;
  
  arma::mat alpha = arma::reshape(ag.rows(0, k_a - 1), k, r);
  
  arma::mat g = arma::zeros<arma::mat>(1, 1);
  if (incl_x) {
    g = arma::reshape(ag.rows(k_a, k_ag - 1), k, k_x);
    y = y - g * x;
  }
  
  arma::mat A = alpha * arma::inv(arma::sqrtmat_sympd(arma::trans(alpha) * alpha));
  arma::mat S_B_post = arma::kron(arma::trans(A) * Sigma_i * A, ect * arma::trans(ect));
  S_B_post = S_B_post + arma::kron(arma::trans(A) * G_i * A, v_i * P_tau_i);
  S_B_post = arma::inv(S_B_post);
  arma::vec mu_B_post = S_B_post * arma::reshape(ect * arma::trans(y) * Sigma_i * A, k_b, 1);
  
  arma::mat U_B;
  arma::mat V_B;
  arma::vec s_B;
  arma::svd(U_B, s_B, V_B, S_B_post);
  arma::mat B_sqrt = U_B * arma::diagmat(sqrt(s_B)) * arma::trans(V_B);
  arma::vec z_B = arma::randn<arma::vec>(k_b);
  arma::mat B = arma::reshape(mu_B_post + B_sqrt * z_B, ect.n_rows, r);
  
  arma::mat BB_sqrt = arma::sqrtmat_sympd(arma::trans(B) * B);
  alpha = A * BB_sqrt;
  beta = B * arma::inv(BB_sqrt);
  arma::mat Pi = alpha * arma::trans(beta);
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("Gamma") = g);
}