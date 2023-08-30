#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::interfaces(r, cpp)]]

//' Durbin and Koopman Simulation Smoother
//' 
//' An implementation of the Kalman filter and backward smoothing
//' algorithm proposed by Durbin and Koopman (2002).
//' 
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param sigma_u the constant \eqn{K \times K} error variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{KT \times K} can be specified.
//' @param sigma_v the constant \eqn{M \times M} coefficient variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{MT \times M} can be specified.
//' @param B an \eqn{M \times M} autocorrelation matrix of the transition equation.
//' @param a_init an M-dimensional vector of initial states.
//' @param P_init an \eqn{M \times M} variance-covariance matrix of the initial states.
//' 
//' @details The function uses algorithm 2 from Durbin and Koopman (2002) to produce
//' a draw of the state vector \eqn{a_t} for \eqn{t = 1,...,T} for a state space model
//' with measurement equation
//' \deqn{y_t = Z_t a_t + u_t}
//' and transition equation 
//' \deqn{a_{t + 1} = B_t a_{t} + v_t,}
//' where \eqn{u_t \sim N(0, \Sigma_{u,t})} and \eqn{v_t \sim N(0, \Sigma_{v,t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = z_t^{\prime} \otimes I_K} is a \eqn{K \times M} matrix of regressors with
//' \eqn{z_t} as a vector of regressors.
//' 
//' The algorithm takes into account Jarociński (2015), where a possible missunderstanding
//' in the implementation of the algorithm of Durbin and Koopman (2002) is pointed out. Following
//' that note the function sets the mean of the initial state to zero in the first step of the algorithm.
//' 
//' @return A \eqn{M \times T+1} matrix of state vector draws.
//' 
//' @examples
//' 
//' # Load data
//' data("e1")
//' data <- diff(log(e1))
//' 
//' # Generate model data
//' temp <- gen_var(data, p = 2, deterministic = "const")
//' y <- t(temp$data$Y)
//' z <- temp$data$SUR
//' k <- nrow(y)
//' tt <- ncol(y)
//' m <- ncol(z)
//' 
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(0.1, m)
//' 
//' a_Q <- diag(.0001, m)
//' 
//' # Initial value of Sigma
//' sigma <- tcrossprod(y) / tt
//' sigma_i <- solve(sigma)
//' 
//' # Initial values for Kalman filter
//' y_init <- y * 0
//' a_filter <- matrix(0, m, tt + 1)
//' 
//' # Initialise the Kalman filter
//' for (i in 1:tt) {
//'   y_init[, i] <- y[, i] - z[(i - 1) * k + 1:k,] %*% a_filter[, i]
//' }
//' a_init <- post_normal_sur(y = y_init, z = z, sigma_i = sigma_i,
//'                           a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//' y_filter <- matrix(y) - z %*% a_init
//' y_filter <- matrix(y_filter, k) # Reshape
//' 
//' # Kalman filter and backward smoother
//' a_filter <- kalman_dk(y = y_filter, z = z, sigma_u = sigma,
//'                       sigma_v = a_Q, B = diag(1, m),
//'                       a_init = matrix(0, m), P_init = a_Q)
//'                       
//' a <- a_filter + matrix(a_init, m, tt + 1)
//' 
//' @references
//' 
//' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for
//' state space time series analysis. \emph{Biometrika, 89}(3), 603--615.
//' 
//' Jarociński, M. (2015). A note on implementing the Durbin and Koopman simulation
//' smoother. \emph{Computational Statistics and Data Analysis, 91}, 1--3.
//' \doi{10.1016/j.csda.2015.05.001}
//' 
// [[Rcpp::export]]
arma::mat kalman_dk(arma::mat y, arma::mat z,
                    arma::mat sigma_u, arma::mat sigma_v,
                    arma::mat B, arma::vec a_init, arma::mat P_init) {
  
  arma::uword k = y.n_rows;
  int t = y.n_cols;
  arma::uword nvars = z.n_cols;
  
  arma::mat sigma_u_temp = arma::zeros<arma::mat>(k * t, k);
  if (sigma_u.n_rows == k){
    for (int i = 0; i < t; i++){
      sigma_u_temp.rows(i * k, (i + 1) * k - 1) = sigma_u;
    }
    sigma_u = sigma_u_temp;
  }
  
  arma::mat sigma_v_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (sigma_v.n_rows == nvars){
    for (int i = 0; i < t; i++){
      sigma_v_temp.rows(i * nvars, (i + 1) * nvars - 1) = sigma_v;
    }
    sigma_v = sigma_v_temp;
  }
  
  arma::mat B_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (B.n_rows == nvars){
    for (int i = 0; i < t; i++){
      B_temp.rows(i * nvars, (i + 1) * nvars - 1) = B;
    }
    B = B_temp;
  }
  
  arma::mat yplus = y * 0;
  arma::mat aplus = arma::zeros<arma::mat>(nvars, t + 1);
  
  arma::mat U;
  arma::vec s;
  
  // Step 1
  arma::eig_sym(s, U, P_init);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
  aplus.col(0) = A * arma::randn<arma::vec>(nvars); // cf. Jarocinski (2015)
  
  int p1, p2, pA1, pA2;
  for (int i = 0; i < t; i++){
    p1 = i * k;
    p2 = (i + 1) * k - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    
    arma::eig_sym(s, U, sigma_u.rows(p1, p2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
    yplus.col(i) = z.rows(p1, p2) * aplus.col(i) + A * arma::randn<arma::vec>(k);
    
    arma::eig_sym(s, U, sigma_v.rows(pA1, pA2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
    aplus.col(i + 1) = B.rows(pA1, pA2) * aplus.col(i) + A * arma::randn<arma::vec>(nvars);
  }
  
  // Step 2
  arma::mat ystar = y - yplus;
  arma::mat a = arma::zeros<arma::mat>(nvars, t + 1);
  a.col(0) = a_init; // cf. DK (2002, p.606)
  arma::mat P = P_init;
  arma::mat v = y * 0;
  arma::mat Fi = arma::zeros<arma::mat>(k * t, k);
  arma::mat K = arma::zeros<arma::mat>(nvars * t, k);
  arma::mat L = arma::zeros<arma::mat>(nvars * t, nvars);
  for (int i = 0; i < t ; i++){
    p1 = i * k;
    p2 = (i + 1) * k - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    v.col(i) = ystar.col(i) - z.rows(p1, p2) * a.col(i);
    Fi.rows(p1, p2) = arma::inv(z.rows(p1, p2) * P * arma::trans(z.rows(p1, p2)) + sigma_u.rows(p1, p2));
    K.rows(pA1, pA2) = B.rows(pA1, pA2) * P * arma::trans(z.rows(p1, p2)) * Fi.rows(p1, p2);
    L.rows(pA1, pA2) = B.rows(pA1, pA2) - K.rows(pA1, pA2) * z.rows(p1, p2);
    a.col(i + 1) = B.rows(pA1, pA2) * a.col(i) + K.rows(pA1, pA2) * v.col(i);
    P = B.rows(pA1, pA2) * P * arma::trans(L.rows(pA1, pA2)) + sigma_v.rows(pA1, pA2);
  }
  
  arma::mat r = arma::zeros<arma::mat>(nvars, t);
  for (int i = (t - 1); i > 0; i--){
    r.col(i - 1) = arma::trans(z.rows(i * k, (i + 1) * k - 1)) * Fi.rows(i * k, (i + 1) * k - 1) * v.col(i) + arma::trans(L.rows(i * nvars, (i + 1) * nvars - 1)) * r.col(i);
  }
  arma::vec r0 = arma::trans(z.rows(0, k - 1)) * Fi.rows(0, k - 1) * v.col(0) + arma::trans(L.rows(0, nvars - 1)) * r.col(0);
  
  a.col(0) = a_init + P_init * r0;
  for (int i = 0; i < t; i++){
    a.col(i + 1) = B.rows(i * nvars, (i + 1) * nvars - 1) * a.col(i) + sigma_v.rows(i * nvars, (i + 1) * nvars - 1) * r.col(i);
  }
  
  // Step 3
  return a + aplus;
}