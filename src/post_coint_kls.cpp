#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw for Cointegration Models
//' 
//' Produces a draw of coefficients for cointegration models with a prior on
//' the cointegration space as proposed in Koop et al. (2010) and a draw of
//' non-cointegration coefficients from a normal density.
//' 
//' @param y a \eqn{K \times T} matrix of differenced endogenous variables.
//' @param beta a \eqn{M \times r} cointegration matrix \eqn{\beta}.
//' @param w a \eqn{M \times T} matrix of variables in the cointegration term.
//' @param x  a \eqn{N \times T} matrix of differenced regressors and unrestricted deterministic terms.
//' @param sigma_i an inverse of the \eqn{K \times K} variance-covariance matrix.
//' @param v_i a numeric between 0 and 1 specifying the shrinkage of the cointegration space prior.
//' @param p_tau_i an inverted \eqn{M \times M} matrix specifying the central location
//' of the cointegration space prior of \eqn{sp(\beta)}.
//' @param g_i a \eqn{K \times K} matrix.
//' @param gamma_mu_prior a \eqn{KN \times 1} prior mean vector of non-cointegration coefficients.
//' @param gamma_v_i_prior an inverted \eqn{KN \times KN} prior covariance matrix of non-cointegration coefficients.
//' 
//' @details The function produces posterior draws of the coefficient
//' matrices \eqn{\alpha}, \eqn{\beta} and \eqn{\Gamma} for the model
//' \deqn{y_{t} = \alpha \beta^{\prime} w_{t-1} + \Gamma z_{t} + u_{t},}
//' where \eqn{y_{t}} is a K-dimensional vector of differenced endogenous variables.
//' \eqn{w_{t}} is an \eqn{M \times 1} vector of variables in the cointegration term,
//' which include lagged values of endogenous and exogenous variables in levels and
//' restricted deterministic terms. \eqn{z_{t}} is an N-dimensional vector of
//' differenced endogenous and exogenous explanatory variabes as well as unrestricted
//' deterministic terms. The error term is \eqn{u_t \sim \Sigma}.
//' 
//' Draws of the loading matrix \eqn{\alpha} are obtained using the prior on the cointegration space
//' as proposed in Koop et al. (2010). The posterior covariance matrix is
//' \deqn{\overline{V}_{\alpha} = \left[\left(v^{-1} (\beta^{\prime} P_{\tau}^{-1} \beta) \otimes G_{-1}\right) + \left(ZZ^{\prime} \otimes \Sigma^{-1} \right) \right]^{-1}}
//' and the posterior mean by
//' \deqn{\overline{\alpha} = \overline{V}_{\alpha} + vec(\Sigma^{-1} Y Z^{\prime}),}
//' where \eqn{Y} is a \eqn{K \times T} matrix of differenced endogenous variables and
//' \eqn{Z = \beta^{\prime} W} with \eqn{W} as an \eqn{M \times T} matrix of
//' variables in the cointegration term.
//' 
//' For a given prior mean vector \eqn{\underline{\Gamma}} and prior covariance matrix \eqn{\underline{V_{\Gamma}}}
//' the posterior covariance matrix of non-cointegration coefficients in \eqn{\Gamma} is obtained by
//' \deqn{\overline{V}_{\Gamma} = \left[ \underline{V}_{\Gamma}^{-1} + \left(X X^{\prime} \otimes \Sigma^{-1} \right) \right]^{-1}}
//' and the posterior mean by
//' \deqn{\overline{\Gamma} = \overline{V}_{\Gamma} \left[ \underline{V}_{\Gamma}^{-1} \underline{\Gamma} + vec(\Sigma^{-1} Y X^{\prime}) \right],}
//' where \eqn{X} is an \eqn{M \times T} matrix of
//' explanatory variables, which do not enter the cointegration term.
//' 
//' Draws of the cointegration matrix \eqn{\beta} are obtained using the prior on the cointegration space
//' as proposed in Koop et al. (2010). The posterior covariance matrix of the unrestricted cointegration
//' matrix \eqn{B} is
//' \deqn{\overline{V}_{B} = \left[\left(A^{\prime} G^{-1} A \otimes v^{-1} P_{\tau}^{-1} \right) + \left(A^{\prime} \Sigma^{-1} A \otimes WW^{\prime} \right) \right]^{-1}}
//' and the posterior mean by
//' \deqn{\overline{B} = \overline{V}_{B} + vec(W Y_{B}^{-1} \Sigma^{-1} A),}
//' where \eqn{Y_{B} = Y - \Gamma X} and \eqn{A = \alpha (\alpha^{\prime} \alpha)^{-\frac{1}{2}}}.
//' 
//' The final draws of \eqn{\alpha} and \eqn{\beta} are calculated using
//' \eqn{\beta = B (B^{\prime} B)^{-\frac{1}{2}}} and
//' \eqn{\alpha = A (B^{\prime} B)^{\frac{1}{2}}}.
//' 
//' @return A named list containing the following elements:
//' \item{alpha}{a draw of the \eqn{K \times r} loading matrix.}
//' \item{beta}{a draw of the \eqn{M \times r} cointegration matrix.}
//' \item{Pi}{a draw of the \eqn{K \times M} cointegration matrix \eqn{\Pi = \alpha \beta^{\prime}}.}
//' \item{Gamma}{a draw of the \eqn{K \times N} coefficient matrix for non-cointegration parameters.}
//' 
//' @examples
//' # Prepare data
//' data("e6")
//' temp <- gen_vec(e6, p = 1)
//' y <- temp$Y
//' ect <- temp$W
//' 
//' k <- nrow(y)
//' t <- ncol(y)
//' 
//' # Initial value of Sigma
//' sigma <- tcrossprod(y) / t
//' sigma_i <- solve(sigma)
//' 
//' # Initial values of beta
//' beta <- matrix(c(1, -4), k)
//' 
//' # Draw parameters
//' coint <- post_coint_kls(y = y, beta = beta, w = ect, sigma_i = sigma_i,
//'                         v_i = 0, p_tau_i = diag(1, k), g_i = sigma_i)
//' 
//' @references
//' 
//' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
//' simulation for cointegrated models with priors on the cointegration space.
//' \emph{Econometric Reviews, 29}(2), 224-242. \url{https://doi.org/10.1080/07474930903382208}
//' 
// [[Rcpp::export]]
Rcpp::List post_coint_kls(arma::mat y, arma::mat beta, arma::mat w, arma::mat sigma_i,
                          double v_i,
                          arma::mat p_tau_i,
                          arma::mat g_i,
                          Rcpp::Nullable<Rcpp::NumericMatrix> x = R_NilValue, 
                          Rcpp::Nullable<Rcpp::NumericVector> gamma_mu_prior = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericMatrix> gamma_v_i_prior = R_NilValue){

  int k = y.n_rows;
  int r = beta.n_cols;
  int k_a = k * r;
  int k_b = w.n_rows * r;
  
  if (p_tau_i.n_cols != w.n_rows) {
    Rcpp::stop("'p_tau_i' must have the same number of rows and columns as the number of variables in the error correction term.");
  }
  
  int k_x = 0;
  int k_g = 0;
  int k_mu, k_v;
  bool incl_x = false;
  arma::mat X, Gamma_V_i_prior;
  arma::vec Gamma_mu_prior;
  if (x.isNotNull()) {
    incl_x = true;
    X = Rcpp::as<arma::mat>(x);
    Gamma_mu_prior = Rcpp::as<arma::vec>(gamma_mu_prior);
    Gamma_V_i_prior = Rcpp::as<arma::mat>(gamma_v_i_prior);
    k_x = X.n_rows;
    k_g = k * k_x;
    k_mu = Gamma_mu_prior.n_rows;
    k_v = Gamma_V_i_prior.n_rows;
    if (k_g != k_mu) {
      Rcpp::stop("Argument 'gamma_mu_prior' does not contain the required amount of elements.");
    }
    if (k_g != k_v) {
      Rcpp::stop("Argument 'gamma_v_i_prior' does not contain the required amount of elements.");
    }
  }
  
  int k_ag = k_a + k_g;
  arma::mat Z = arma::zeros<arma::mat>(r + k_x, y.n_cols);
  Z.rows(0, r - 1) = arma::trans(beta) * w;
  if (incl_x) {
    Z.rows(r, r + k_x - 1) = X; 
  }
  arma::vec mu_ag_prior = arma::zeros<arma::vec>(k_ag);
  arma::mat V_ag_post = arma::kron(Z * arma::trans(Z), sigma_i);
  arma::mat V_ag_prior = V_ag_post * 0;
  V_ag_prior.submat(0, 0, k_a - 1, k_a - 1) = arma::kron(v_i * (arma::trans(beta) * p_tau_i * beta), g_i);
  if (incl_x) {
    mu_ag_prior.subvec(k_a, k_ag - 1) = Gamma_mu_prior;
    V_ag_prior.submat(k_a, k_a, k_ag - 1, k_ag - 1) = Gamma_V_i_prior; 
  }
  V_ag_post = arma::inv(V_ag_post + V_ag_prior);
  arma::vec mu_ag_post = V_ag_post * (mu_ag_prior + arma::reshape(sigma_i * y * arma::trans(Z), k_ag, 1));

  arma::mat U_ag;
  arma::vec s_ag;
  arma::eig_sym(s_ag, U_ag, V_ag_post);
  arma::mat ag_sqrt = U_ag * arma::diagmat(sqrt(s_ag)) * arma::trans(U_ag);
  arma::vec z_ag = arma::randn<arma::vec>(k_ag);
  arma::mat ag = mu_ag_post + ag_sqrt * z_ag;

  arma::mat alpha = arma::reshape(ag.rows(0, k_a - 1), k, r);

  arma::mat g;
  if (incl_x) {
    g = ag.rows(k_a, k_ag - 1);
    y = y - arma::reshape(g, k, k_x) * X;
  }
  
  // Obtain beta
  arma::mat A = alpha * arma::sqrtmat_sympd(arma::inv(arma::trans(alpha) * alpha));
  arma::mat S_B_post = arma::kron(arma::trans(A) * sigma_i * A, w * arma::trans(w));
  S_B_post = S_B_post + arma::kron(arma::trans(A) * g_i * A, v_i * p_tau_i);
  S_B_post = arma::inv(S_B_post);
  arma::vec mu_B_post = S_B_post * arma::reshape(w * arma::trans(y) * sigma_i * A, k_b, 1);
  
  arma::mat U_B, V_B;
  arma::vec s_B;
  arma::svd(U_B, s_B, V_B, S_B_post);
  arma::mat B_sqrt = U_B * arma::diagmat(sqrt(s_B)) * arma::trans(V_B);
  arma::vec z_B = arma::randn<arma::vec>(k_b);
  arma::mat B = arma::reshape(mu_B_post + B_sqrt * z_B, w.n_rows, r);

  arma::mat BB_sqrt = arma::sqrtmat_sympd(arma::trans(B) * B);
  alpha = A * BB_sqrt;
  beta = B * arma::inv(BB_sqrt);
  arma::mat Pi = alpha * arma::trans(beta);

  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("Gamma") = g);
}