#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Posterior Data Preparation
//' 
//' Generates the input matrix for the simulation of the cointegration matrix \eqn{\beta}.
//' 
//' @param w a \eqn{T \times M} matrix.
//' @param alpha a \eqn{K \times r} matrix of the loading matrix \eqn{\alpha}.
//' @param k integer of the number of endogenous variables in the model.
//' @param r integer with the rank of the cointegration matrix \eqn{\Pi}.
//' @param reparameterise logical indicating if the reparameterisation described
//' in Koop et al. (2010) should be performed.
//' @param tvp logical indicating if the SUR matrix with the values of regressors
//' should be prepared for the estimation of constant or time varying parameters.
//' 
//' @return A sparse block matrix.
//' 
//' @examples
//' 
//' data("e6") # Load data
//' 
//' # Generate model input
//' mod <- create_bvecmodel(e6, p = 4, r = 1,
//'                         const = "unrestricted", seasonal = "unrestricted")
//' 
//' # Obtain input data
//' alpha <- matrix(c(-0.1, 0.16, -0.04, -0.02), 2)
//' w <- mod$data$train$w
//' 
//' # Constant coefficients
//' coint_prepare_sur_data(w, alpha, 2, 1, TRUE, FALSE)
//' 
//' # Time varying coefficients
//' coint_prepare_sur_data(w, alpha, 2, 1, FALSE, TRUE)
//' 
//' 
// [[Rcpp::export(coint_prepare_sur_data)]]
Rcpp::List coint_prepare_sur_data(const arma::mat w, arma::mat alpha, const int k, const int r, const bool reparameterise, const bool tvp) {
  
  const int tt = w.n_rows;
  const int n_beta = r * w.n_cols;
  
  // Reparameterise alpha
  alpha = arma::reshape(alpha, k, r);
  arma::mat alpha_new;
  if (reparameterise) {
    alpha_new = alpha * arma::solve(arma::sqrtmat_sympd(alpha.t() * alpha), arma::eye(r, r)); 
  } else {
    alpha_new = alpha;
  }
  
  // Generate SUR data for simulation of beta
  arma::sp_mat z;
  
  if (tvp) {
    z = arma::zeros<arma::sp_mat>(k * tt, n_beta * tt);
    for (int i = 0; i < tt; i++){
      z.submat(i * k, i * n_beta, (i + 1) * k - 1, (i + 1) * n_beta - 1) = arma::kron(alpha_new, w.row(i));
    } 
  } else {
    z = arma::zeros<arma::sp_mat>(k * tt, n_beta);
    for (int i = 0; i < tt; i++){
      z.rows(i * k, (i + 1) * k - 1) = arma::kron(alpha_new, w.row(i));
    } 
  }
  
  return Rcpp::List::create(Rcpp::Named("z") = z,
                            Rcpp::Named("alpha") = alpha_new);
}

/*** R

data("e6") # Load data

# Generate model input
mod <- create_bvecmodel(e6, p = 4, r = 1,
                        const = "unrestricted", seasonal = "unrestricted")

# Obtain input data
alpha <- matrix(c(-0.1, 0.16, -0.04, -0.02), 2)
w <- mod$data$train$w

# Constant coefficients
coint_prepare_sur_data(w, alpha, 2, 1, TRUE, FALSE)

# Time varying coefficients
coint_prepare_sur_data(w, alpha, 2, 1, FALSE, TRUE)
*/