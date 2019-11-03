#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw from a Normal Distribution
//' 
//' Produces a draw of coefficients from a normal posterior density for
//' a model with seemingly unrelated regresssions (SUR).
//' 
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param sigma_i the inverse of the constant \eqn{K \times K} error variance-covariance matrix.
//' For time varying variance-covariance matrics a \eqn{KT \times K} can be provided.
//' @param a_prior a \eqn{M x 1} numeric vector of prior means.
//' @param v_i_prior the inverse of the \eqn{M x M} prior covariance matrix.
//' 
//' @details The function produces a posterior draw of the coefficient vector \eqn{a} for the model
//' \deqn{y_{t} = Z_{t} a + u_{t},}
//' where \eqn{u_t \sim N(0, \Sigma_{t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = z_t^{\prime} \otimes I_K} is a \eqn{K \times KM} matrix of regressors with
//' \eqn{z_t} as a vector of regressors.
//' 
//' For a given prior mean vector \eqn{\underline{a}} and prior covariance matrix \eqn{\underline{V}}
//' the posterior covariance matrix is obtained by
//' \deqn{\overline{V} = \left[ \underline{V}^{-1} + \sum_{t=1}^{T} Z_{t}^{\prime} \Sigma_{t}^{-1} Z_{t} \right]^{-1}}
//' and the posterior mean by
//' \deqn{\overline{a} = \overline{V} \left[ \underline{V}^{-1} \underline{a} + \sum_{t=1}^{T} Z_{t}^{\prime} \Sigma_{t}^{-1} y_{t}  \right].}
//' 
//' @examples
//' # Prepare data
//' data("e1")
//' data <- diff(log(e1))
//' temp <- gen_var(data, p = 2, deterministic = "const")
//' y <- temp$Y
//' x <- temp$Z
//' k <- nrow(y)
//' z <- kronecker(t(x), diag(1, k))
//' t <- ncol(y)
//' m <- k * nrow(x)
//' 
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(0.1, m)
//' 
//' # Initial value of inverse Sigma
//' sigma_i <- solve(tcrossprod(y) / t)
//'
//' # Draw parameters
//' a <- post_normal_sur(y = y, z = z, sigma_i = sigma_i,
//'                      a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//' 
//' @return A vector.
//' 
// [[Rcpp::export]]
arma::vec post_normal_sur(arma::mat y ,arma::mat z, arma::mat sigma_i,
                          arma::vec a_prior, arma::mat v_i_prior) {
  
  arma::uword n = y.n_rows;
  int t = y.n_cols;
  int nvars = z.n_cols;
  int n_mu = a_prior.n_rows;
  int n_v = v_i_prior.n_rows;
  if (nvars != n_mu) {
    Rcpp::stop("Argument 'a_prior' does not contain the required amount of elements.");
  }
  if (nvars != n_v) {
    Rcpp::stop("Argument 'v_i_prior' does not contain the required amount of elements.");
  }
  
  bool const_var = true;
  if (sigma_i.n_rows > n) {
    const_var = false;
  }
  
  arma::mat ZHi = arma::zeros<arma::mat>(nvars, n * t);
  arma::mat ZHZ = arma::zeros<arma::mat>(nvars, nvars);
  arma::mat ZHy = arma::zeros<arma::mat>(nvars, 1);
  
  arma::mat S_i;
  if (const_var) {
    ZHi = arma::trans(z) * arma::kron(arma::eye<arma::mat>(t, t), sigma_i);
  } else {
    S_i = arma::zeros<arma::mat>(n * t, n * t);
    for (int i = 0; i < t; i++){
      S_i.submat(i * n, i * n, (i + 1) * n - 1, (i + 1) * n - 1) = sigma_i.rows(i * n, (i + 1) * n - 1);
    }
    ZHi = arma::trans(z) * S_i;
  }
  ZHZ = ZHi * z;
  ZHy = ZHi * arma::vectorise(y);
  
  arma::mat V_post = arma::inv(v_i_prior + ZHZ);
  arma::vec mu_post = V_post * (v_i_prior * a_prior + ZHy);

  arma::vec s;
  arma::mat U;
  arma::eig_sym(s, U, V_post);
  
  return mu_post + U * arma::diagmat(sqrt(s)) * arma::trans(U) * arma::randn<arma::vec>(nvars);
}