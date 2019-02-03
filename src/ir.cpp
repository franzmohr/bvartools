#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Impulse Response
//' 
//' Produces different types of impulse responses.
//' 
//' @param A \eqn{K \times Kp} matrix of coefficients, where \eqn{K} is the number of endogenous
//' variables and \eqn{p} is the number of lags.
//' @param h integer specifying the steps.
//' @param type character specifying the type of the impulse response.
//' @param impulse numeric specifying the position of the impulse variable.
//' @param response numeric specifying the position of the response variable.
//' 
//' @details The function produces different types of impulse responses for the VAR model
//' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
//' with \eqn{u_t \sim N(0, \Sigma)}.
//' 
//' Forecast error impulse responses \eqn{\Phi_i} are obtained by recursions
//' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i-j} A_j,   i = 1, 2,...,}
//' with \eqn{\Phi_0 = I_K}.
//' 
//' Orthogonalised impulse responses \eqn{\Theta^o_i} are calculated as \eqn{\Theta^o_i = \Phi_i P,}
//' where P is the lower triangular Choleski decomposition of \eqn{\Sigma}.
//' 
//' Structural impulse responses \eqn{\Theta^s_i} are calculated as \eqn{\Theta^s_i = \Phi_i A_0^{-1}}.
//' 
//' Generalised impulse responses for variable \eqn{j}, i.e. \eqn{\Theta^g_ji} are calculated as
//' \eqn{\Theta^g_ji = \sigma_{jj}^{-1/2} \Phi_i \Sigma e_j}, where \eqn{\sigma_{jj}} is the variance
//' of the \eqn{j^{th}} diagonal element of \eqn{\Sigma} and \eqn{e_i} is a selection vector containing
//' one in its \eqn{j^{th}} element and zero otherwise.
//' 
//' @return A \eqn{K(h + 1) \times K} matrix.
//' 
//' @references
//' 
//' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
//' 
// [[Rcpp::export]]
arma::vec ir(Rcpp::List A, int h, std::string type, int impulse, int response) {
  arma::mat coef = Rcpp::as<arma::mat>(A["A"]);
  
  int k = coef.n_rows;
  int p = coef.n_cols / k;
  
  std::string fe ("feir");
  std::string oir ("oir");
  std::string sir ("sir");
  std::string gir ("gir");
  
  if (h < p) {
    p = h * k;
  } else {
    p = coef.n_cols;
  }
  
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = coef.cols(0, p - 1); 
  
  arma::mat A0, P, Sigma, temp;
  
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp;
  arma::vec theta = arma::zeros<arma::vec>(h + 1);
  
  if (type == fe) {
    theta(0) = arma::as_scalar(phi(response - 1, impulse - 1));
  }
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"])));
  }
  if (type == sir) {
    A0 = Rcpp::as<arma::mat>(A["A0"]);
    P = arma::inv(A0);
  }
  if (type == gir) {
    Sigma = Rcpp::as<arma::mat>(A["Sigma"]);
    P = Sigma / sqrt(arma::as_scalar(Sigma(impulse - 1, impulse - 1)));
  }
  if (type != fe) {
    temp = phi_temp * P;
    theta(0) = arma::as_scalar(temp(response - 1, impulse - 1));
  }
  
  for (int i = 1; i <= h; i++) {
    phi_temp.zeros();
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp;
    
    if (type == fe) {
      theta(i) = arma::as_scalar(phi_temp(response - 1, impulse - 1));
    } else {
      temp = phi_temp * P;
      theta(i) = arma::as_scalar(temp(response - 1, impulse - 1));
    }
  }

  return theta;
}