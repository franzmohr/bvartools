#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Bayesian Variable Selection
//' 
//' \code{bvs} employs Bayesian variable selection as proposed by Korobilis (2013)
//' to produce a vector of inclusion parameters for the coefficient matrix
//' of a VAR model.
//' 
//' @param y a \eqn{K \times T} matrix of the endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param a an M-dimensional vector of parameter draws. If time varying parameters are used,
//' an \eqn{M \times T} coefficient matrix can be provided.
//' @param lambda an \eqn{M \times M} inclusion matrix that should be updated.
//' @param sigma_i the inverse variance-covariance matrix. If the variance-covariance matrix
//' is time varying, a \eqn{KT \times K} matrix can be provided.
//' @param prob_prior an M-dimensional vector of prior inclusion probabilities.
//' @param include an integer vector specifying the positions of variables, which should be
//' included in the BVS algorithm. If \code{NULL} (default), BVS will be applied to all variables.
//' 
//' @details The function employs Bayesian variable selection as proposed
//' by Korobilis (2013) to produce a vector of inclusion parameters, which are
//' the diagonal elements of the inclusion matrix \eqn{\Lambda} for the VAR model
//' \deqn{y_t = Z_t \Lambda a_t + u_t,}
//' where \eqn{u_t \sim N(0, \Sigma_{t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = x_t^{\prime} \otimes I_K} is a \eqn{K \times M} matrix of regressors with
//' \eqn{x_t} as a vector of regressors.
//' 
//' @return A matrix of inclusion parameters on its diagonal.
//' 
//' @examples
//' 
//' # Load data
//' data("e1")
//' data <- diff(log(e1)) * 100
//' 
//' # Generate model data
//' temp <- gen_var(data, p = 2, deterministic = "const")
//' 
//' y <- t(temp$data$Y)
//' z <- temp$data$SUR
//' 
//' tt <- ncol(y)
//' m <- ncol(z)
//' 
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(0.1, m)
//' 
//' # Prior for inclusion parameter
//' prob_prior <- matrix(0.5, m)
//' 
//' # Initial value of Sigma
//' sigma <- tcrossprod(y) / tt
//' sigma_i <- solve(sigma)
//' 
//' lambda <- diag(1, m)
//' 
//' z_bvs <- z %*% lambda
//' 
//' a <- post_normal_sur(y = y, z = z_bvs, sigma_i = sigma_i,
//'                      a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//'
//' lambda <- bvs(y = y, z = z, a = a, lambda = lambda,
//'               sigma_i = sigma_i, prob_prior = prob_prior)
//' 
//' @references
//' 
//' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection. \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
//' 
// [[Rcpp::export]]
arma::mat bvs(arma::mat y, arma::mat z, arma::mat a, arma::mat lambda, arma::mat sigma_i,
              arma::vec prob_prior, Rcpp::Nullable<Rcpp::IntegerVector> include = R_NilValue)
{
  arma::vec lpr_prior_0 = arma::log(1 - prob_prior);
  arma::vec lpr_prior_1 = arma::log(prob_prior);
  int k = y.n_rows;
  int m = a.n_rows;
  int t = y.n_cols;
  int nomega = sigma_i.n_rows;
  
  bool const_par = true;
  if (a.n_cols > 1){
    const_par = false;
  }
  
  bool sv = false;
  if (nomega > k) {
    sv = true;
  }
  
  arma::mat S_i = arma::zeros<arma::mat>(k * t, k * t);
  if (sv) {
    for (int i = 0; i < t; i++) {
      S_i.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1 ) = sigma_i.rows(i * k, (i + 1) * k - 1);
    }
  } else {
    S_i = kron(arma::eye<arma::mat>(t, t), sigma_i);  
  }
  
  arma::vec yvec = vectorise(y);
  arma::vec g = arma::zeros<arma::vec>(1);
  arma::mat AG = lambda * a;
  arma::mat theta0 = AG;
  arma::mat theta1 = AG;
  arma::vec l0_res = arma::zeros<arma::vec>(k * t);
  arma::vec l1_res = arma::zeros<arma::vec>(k * t);
  
  double l0 = 0;
  double l1 = 0;
  double bayes = 0;
  double bayes_rand = 0;
  
  arma::vec pos_res(m);
  for (int l = 0; l < m; l++) {
    pos_res(l) = l;
  }
  arma::uvec ex;
  if (include.isNotNull()) {
    ex = Rcpp::as<arma::uvec>(include) - 1;
    pos_res = pos_res.elem(ex);
  }
  
  pos_res = shuffle(pos_res);
  int pos_size = size(pos_res)(0);
  int var = 0;
  for (int j = 0; j < pos_size; j++){
    var = pos_res(j);
    g = log(arma::randu<arma::vec>(1));
    if (lambda(var, var) == 1 && g(0) >= lpr_prior_1(var)){continue;}
    if (lambda(var, var) == 0 && g(0) >= lpr_prior_0(var)){continue;}
    if ((lambda(var, var) == 1 && g(0) < lpr_prior_1(var)) || (lambda(var, var) == 0 && g(0) < lpr_prior_0(var))){
      theta0 = AG;
      theta1 = AG;
      theta1.row(var) = a.row(var);
      if (const_par) {
        theta0.row(var) = 0;
        l0_res = yvec - z * theta0;
        l1_res = yvec - z * theta1;
      } else {
        theta0.row(var) = arma::zeros<arma::mat>(1, t);
        for (int i = 0; i < t; i++){
          l0_res.rows(i * k, (i + 1) * k - 1) = yvec.rows(i * k, (i + 1) * k - 1) - z.rows(i * k, (i + 1) * k - 1) * theta0.col(i);
          l1_res.rows(i * k, (i + 1) * k - 1) = yvec.rows(i * k, (i + 1) * k - 1) - z.rows(i * k, (i + 1) * k - 1) * theta1.col(i);
        }
      }
      l0 = -.5 * arma::as_scalar(trans(l0_res) * S_i * l0_res) + arma::as_scalar(lpr_prior_0(var));
      l1 = -.5 * arma::as_scalar(trans(l1_res) * S_i * l1_res) + arma::as_scalar(lpr_prior_1(var));
      bayes = l1 - l0;
      bayes_rand = log(arma::as_scalar(arma::randu<arma::vec>(1)));
      if (bayes >= bayes_rand){
        lambda(var, var) = 1;
      } else {
        lambda(var, var) = 0;
      }
    }
  }
  return lambda;
}