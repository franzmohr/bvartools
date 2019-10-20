#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

//' Algorithm of Chan and Jeliazkov (2009)
//' 
//' A simplified implementation of the algorithm of Chan and Jeliazkov (2009).
//' 
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param sigma_i the inverse of the constant \eqn{K \times K} error variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{KT \times K} can be specified.
//' @param q_i the inverse of the constant \eqn{M \times M} coefficient variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{MT \times M} can be specified.
//' @param a0 an M-dimensional vector of initial states.
//' 
//' @details The function uses a simplified version of the algorithm of of Chan and Jeliazkov (2009)
//' to produce a draw of the state vector \eqn{a_t} for \eqn{t = 1,...,T} for a state space model
//' with measurement equation
//' \deqn{y_t = Z_t a_t + u_t}
//' and transition equation 
//' \deqn{a_{t} = a_{t - 1} + v_t,}
//' where \eqn{u_t \sim N(0, \Sigma_{t})} and \eqn{v_t \sim N(0, Q_{t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = z_t^{\prime} \otimes I_K} is a \eqn{K \times M} matrix of regressors with
//' \eqn{z_t} as a vector of regressors.
//' 
//' The implementation follows the depiction in chapter 20 of Chan, Koop, Poirier and Tobias (2019).
//' 
//' @return A \eqn{M \times T} matrix of state vector draws.
//' 
//' @examples
//' # Load data and transform it
//' data("e1")
//' e1 <- diff(log(e1[, c("cons", "income")]))
//' 
//' # Create model
//' model_data <- gen_var(e1, p = 2, deterministic = "const")
//' y <- matrix(model_data$Y["cons", ], 1)
//' x <- rbind(model_data$Y["income", ], model_data$Z)
//' 
//' tt <- NCOL(y) # Number of periods
//' k <- NROW(y) # Number of endogenous variables
//' 
//' z <- kronecker(t(x), diag(1, k)) # Matrix of regressors
//' 
//' # Initial values
//' b0 <- matrix(tcrossprod(y, x) %*% solve(tcrossprod(x)))
//' epsilon <- matrix(matrix(y) - z %*% b0, k)
//' sigma <- tcrossprod(epsilon) / tt
//' sigma_i <- solve(sigma)
//' q_i <- diag(100, k * nrow(x))
//' 
//' b <- chan_jeliazkov(y = y, z = z, sigma_i = sigma_i,
//'                     q_i = q_i, a0 = b0)
//' 
//' @references
//' 
//' Chan, J. C. C., & Jeliazkov, I. (2009). Efficient simulation and integrated likelihood
//' estimation in state space models. \emph{International Journal of Mathematical Modelling
//' and Numerical Optimisation, 1}(1/2), 101--120. \url{https://doi.org/10.1504/ijmmno.2009.030090}
//' 
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
// [[Rcpp::export]]
arma::mat chan_jeliazkov(arma::mat y, arma::mat z,
                    arma::mat sigma_i, arma::mat q_i,
                    arma::vec a0) {
  
  arma::uword k = y.n_rows;
  int tt = y.n_cols;
  arma::uword nvars = z.n_cols;
  
  // Prepare Z
  arma::sp_mat Z = arma::zeros<arma::sp_mat>(k * tt, nvars * tt);
  for (int i = 0; i < tt; i++){
    Z.submat(i * k, i * nvars, (i + 1) * k - 1, (i + 1) * nvars - 1) = z.rows(i * k, (i + 1) * k - 1);
  }
  
  // Prepare Sigma
  arma::sp_mat S_i;
  if (sigma_i.n_rows == k){
    S_i = arma::sp_mat(arma::kron(arma::eye(tt, tt), sigma_i));
  } else {
    S_i = arma::zeros<arma::sp_mat>(k * tt, k * tt);
    for (int i = 0; i < tt; i++){
      S_i.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1) = sigma_i.rows(i * k, (i + 1) * k - 1);
    }
  }
  
  // Prepare Q
  arma::sp_mat Q_i;
  if (q_i.n_rows == nvars){
    Q_i = arma::sp_mat(arma::kron(arma::eye(tt, tt), q_i));
  } else {
    Q_i = arma::zeros<arma::sp_mat>(nvars * tt, nvars * tt);
    for (int i = 0; i < tt; i++){
      Q_i.submat(i * nvars, i * nvars, (i + 1) * nvars - 1, (i + 1) * nvars - 1) = q_i.rows(i * nvars, (i + 1) * nvars - 1);
    }
  }
  
  // Prepare B
  arma::sp_mat B = arma::eye<arma::sp_mat>(tt * nvars, tt * nvars);
  B.diag(-nvars) += -1;
  
  arma::sp_mat HQiH = arma::trans(B) * Q_i * B;
  arma::sp_mat Zsi = arma::trans(Z) * S_i;
  arma::sp_mat K = HQiH + Zsi * Z;
  arma::vec b0 = arma::kron(arma::ones<arma::vec>(tt), a0);
  arma::vec mu = arma::spsolve(K, HQiH * b0 + Zsi * arma::reshape(y, k * tt, 1), "lapack");
  arma::mat Kd(K);
  arma::sp_mat A = arma::sp_mat(arma::trans(arma::chol(Kd, "lower")));
  arma::vec a = mu + arma::spsolve(A, arma::randn<arma::vec>(nvars * tt), "lapack");
  return arma::reshape(a, nvars, tt);
}