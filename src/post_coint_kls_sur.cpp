#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw for Cointegration Models
//' 
//' Produces a draw of coefficients for cointegration models in SUR form with a prior on
//' the cointegration space as proposed in Koop et al. (2010) and a draw of
//' non-cointegration coefficients from a normal density.
//' 
//' @param y a \eqn{K \times T} matrix of differenced endogenous variables.
//' @param beta a \eqn{M \times r} cointegration matrix \eqn{\beta}, where \eqn{\beta^{\prime} \beta = I}.
//' @param w a \eqn{M \times T} matrix of variables in the cointegration term.
//' @param x  a \eqn{KT \times NK} matrix of differenced regressors and unrestricted deterministic terms.
//' @param sigma_i the inverse of the constant \eqn{K \times K} error variance-covariance matrix.
//' For time varying variance-covariance matrics a \eqn{KT \times K} can be provided.
//' @param v_i a numeric between 0 and 1 specifying the shrinkage of the cointegration space prior.
//' @param p_tau_i an inverted \eqn{M \times M} matrix specifying the central location
//' of the cointegration space prior of \eqn{sp(\beta)}.
//' @param g_i a \eqn{K \times K} or \eqn{KT \times K} matrix. If the matrix is \eqn{KT \times K},
//' the function will automatically produce a \eqn{K \times K} matrix containing the means of the
//' time varying \eqn{K \times K} covariance matrix.
//' @param gamma_mu_prior a \eqn{KN \times 1} prior mean vector of non-cointegration coefficients.
//' @param gamma_v_i_prior an inverted \eqn{KN \times KN} prior covariance matrix of
//' non-cointegration coefficients.
//' @param svd logical. If \code{TRUE} the singular value decomposition is used to determine
//' the root of the posterior covariance matrix. Default is \code{FALSE} which means
//' that the eigenvalue decomposition is used.
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
//' 
//' # Load data
//' data("e6")
//' 
//' # Generate model data
//' temp <- gen_vec(e6, p = 1, r = 1)
//' y <- t(temp$data$Y)
//' ect <- t(temp$data$W)
//' 
//' k <- nrow(y) # Endogenous variables
//' tt <- ncol(y) # Number of observations
//' 
//' # Initial value of Sigma
//' sigma <- tcrossprod(y) / tt
//' sigma_i <- solve(sigma)
//' 
//' # Initial values of beta
//' beta <- matrix(c(1, -4), k)
//' 
//' # Draw parameters
//' coint <- post_coint_kls_sur(y = y, beta = beta, w = ect,
//'                             sigma_i = sigma_i, v_i = 0, p_tau_i = diag(1, nrow(ect)),
//'                             g_i = sigma_i)
//' 
//' @references
//' 
//' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
//' simulation for cointegrated models with priors on the cointegration space.
//' \emph{Econometric Reviews, 29}(2), 224-242. \url{https://doi.org/10.1080/07474930903382208}
//' 
// [[Rcpp::export]]
Rcpp::List post_coint_kls_sur(arma::mat y, arma::mat beta, arma::mat w, arma::mat sigma_i,
                              double v_i,
                              arma::mat p_tau_i,
                              arma::mat g_i,
                              Rcpp::Nullable<Rcpp::NumericMatrix> x = R_NilValue, 
                              Rcpp::Nullable<Rcpp::NumericVector> gamma_mu_prior = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericMatrix> gamma_v_i_prior = R_NilValue,
                              bool svd = false){

  if (sigma_i.has_nan()) {
    Rcpp::stop("Argument 'sigma_i' contains NAs.");
  }  
  
  arma::uword k = y.n_rows;
  int t = y.n_cols;
  int r = beta.n_cols;
  int k_a = k * r;
  int k_b = w.n_rows * r;
  
  if (p_tau_i.n_cols != w.n_rows) {
    Rcpp::stop("'p_tau_i' must have the same number of rows and columns as the number of variables in the error correction term.");
  }
  
  bool const_var = true;
  arma::mat S_i = arma::zeros<arma::mat>(k * t, k * t);
  if (sigma_i.n_rows > k) {
    const_var = false;
  }
  
  arma::mat G_i = g_i;
  if (g_i.n_rows > k) {
    G_i = arma::zeros<arma::mat>(k, k);
    for (int i = 0; i < t; i++){
      G_i = G_i + g_i.rows(i * k, (i + 1) * k - 1);
    }
    G_i = G_i / t;
  }
  
  int k_g = 0;
  bool incl_x = false;
  int k_mu, k_v;
  arma::mat X, Gamma_V_i_prior;
  arma::vec Gamma_mu_prior;
  if (x.isNotNull()) {
    incl_x = true;
    X = Rcpp::as<arma::mat>(x);
    Gamma_mu_prior = Rcpp::as<arma::vec>(gamma_mu_prior);
    Gamma_V_i_prior = Rcpp::as<arma::mat>(gamma_v_i_prior);
    k_g = X.n_cols;
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
  
  arma::mat z = arma::zeros<arma::mat>(k * t, k_ag);
  z.cols(0, k_a - 1) = arma::kron(arma::trans(arma::trans(beta) * w), arma::eye<arma::mat>(k, k));
  if (incl_x) {
    z.cols(k_a, k_ag - 1) = X; 
  }
  
  arma::vec mu_ag_prior = arma::zeros<arma::vec>(k_ag);
  arma::mat V_ag_prior = arma::zeros<arma::mat>(k_ag, k_ag);
  V_ag_prior.submat(0, 0, k_a - 1, k_a - 1) = arma::kron(v_i * (arma::trans(beta) * p_tau_i * beta), G_i);
  if (incl_x) {
    mu_ag_prior.subvec(k_a, k_ag - 1) = Gamma_mu_prior;
    V_ag_prior.submat(k_a, k_a, k_ag - 1, k_ag - 1) = Gamma_V_i_prior;
  }
  
  arma::mat ZHi = arma::zeros<arma::mat>(k_ag, k * t);
  arma::mat ZHZ = arma::zeros<arma::mat>(k_ag, k_ag);
  arma::mat ZHy = arma::zeros<arma::mat>(k_ag, 1);
  
  if (const_var) {
    ZHi = arma::trans(z) * arma::kron(arma::eye<arma::mat>(t, t), sigma_i);
  } else {
    for (int i = 0; i < t; i++){
      S_i.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1) = sigma_i.rows(i * k, (i + 1) * k - 1);
    }
    ZHi = arma::trans(z) * S_i;
  }
  ZHZ = ZHi * z;
  ZHy = ZHi * arma::vectorise(y);
  arma::mat V_post_i = arma::inv(V_ag_prior + ZHZ);
  arma::vec mu_post = V_post_i * (V_ag_prior * mu_ag_prior + ZHy);
  
  arma::vec s;
  arma::mat U, V;
  if (svd) {
    arma::svd(U, s, V, V_post_i);
  } else {
    arma::eig_sym(s, U, V_post_i);
    V = U;
  }
  arma::mat ag = mu_post + U * arma::diagmat(sqrt(s)) * arma::trans(V) * arma::randn<arma::vec>(k_ag);
  
  arma::mat alpha = arma::reshape(ag.rows(0, k_a - 1), k, r);
    
  arma::mat g;
  if (incl_x) {
    g = ag.rows(k_a, k_ag - 1);
    y = y - arma::reshape(X * g, k, t);
  }
  
  arma::mat A = alpha * arma::inv(arma::sqrtmat_sympd(arma::trans(alpha) * alpha));
  
  z = arma::zeros<arma::mat>(k * t, k_b);
  for (int i = 0; i < t; i++){
    z.rows(i * k, (i + 1) * k - 1) = arma::kron(A, arma::trans(w.col(i)));
  }
  
  if (const_var) {
    ZHi = arma::trans(z) * arma::kron(arma::eye<arma::mat>(t, t), sigma_i);
  } else {
    ZHi = arma::trans(z) * S_i;
  }
  ZHZ = ZHi * z;
  ZHy = ZHi * arma::vectorise(y);
  
  V_post_i = arma::inv(arma::kron(arma::trans(A) * G_i * A, v_i * p_tau_i) + ZHZ);
  mu_post = V_post_i * ZHy;
  
  arma::svd(U, s, V, V_post_i);
  arma::mat B = mu_post + U * arma::diagmat(arma::sqrt(s)) * arma::trans(V) * arma::randn<arma::vec>(k_b);
  B = arma::reshape(B, w.n_rows, r);
  
  arma::mat BB_sqrt = arma::sqrtmat_sympd(arma::trans(B) * B);
  alpha = A * BB_sqrt;
  beta = B * arma::inv(BB_sqrt);
  arma::mat Pi = alpha * arma::trans(beta);
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("Gamma") = g);
}