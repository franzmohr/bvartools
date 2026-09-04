#include "../inst/include/bvartools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Called as a plain C++ function rather than through this package's own C++
// interface, which would rewind the random number stream. See the header.
#include "sur_const_to_tvp.h"

//' Bayesian Variable Selection
//' 
//' \code{post_bvs} employs Bayesian variable selection as proposed by Korobilis (2013)
//' to produce a vector of inclusion parameters for the coefficient matrix
//' of a VAR model.
//' 
//' @param y a \eqn{KT \times 1} vector of the endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param a an M-dimensional vector of parameter draws. If time varying parameters are used,
//' an \eqn{M \times T} coefficient matrix can be provided.
//' @param k integer of the number of endogenous variables.
//' @param m integer of the number of M 
//' @param lambda an \eqn{M \times M} inclusion matrix that should be updated.
//' @param sigma_i a sparse \eqn{KT \times KT} block diagonal matrix containing
//' the inverse variance-covariance matrix.
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
//' e1 <- diff(log(e1)) * 100
//' 
//' # Generate model input, which uses BVS as variable selection algorithm
//' object <- create_bvarmodel(data = e1, p = 2, deterministic = "const",
//'                            varsel = "bvs")
//' 
//' # Add prior specifications, including the prior inclusion probabilities
//' object <- add_priors(object,
//'                      coef = list(v_i = 1, v_i_det = 1 / 10),
//'                      sigma = list(df = "k", scale = 1),
//'                      varsel = list(inprior = .1))
//' 
//' # Add initial values
//' object <- add_initial_values(object)
//' 
//' # Obtain data, initial values and priors
//' y <- matrix(t(object[["data"]][["train"]][["y"]]))
//' z <- object[["data"]][["train"]][["z"]] # Argument 'z' is taken dense
//' k <- object[["model"]][["k"]]
//' tt <- nrow(object[["data"]][["train"]][["y"]])
//' m <- ncol(z)
//' a <- object[["initial"]][["a"]]
//' prob_prior <- object[["priors"]][["a"]][["inprior"]]
//' 
//' # Arguments 'lambda' and 'sigma_i' have to be sparse
//' lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)
//' 
//' # Initial value of the inverse error covariance matrix
//' u <- matrix(y - z %*% a, k)
//' sigma_i <- Matrix::Matrix(kronecker(diag(1, tt), solve(tcrossprod(u) / tt)),
//'                           sparse = TRUE)
//' 
//' # Draw inclusion parameters
//' post_bvs(y, z, a, k, m, lambda, sigma_i, prob_prior)
//' 
//' @references
//' 
//' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
//' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
//' 
// [[Rcpp::export]]
arma::sp_mat post_bvs(const arma::vec& y, const arma::mat& z, arma::vec& a, arma::uword k, arma::uword m,
                      arma::sp_mat lambda, arma::sp_mat sigma_i, arma::vec prob_prior,
                      Rcpp::Nullable<Rcpp::IntegerVector> include = R_NilValue) {
   
   arma::vec lpr_prior_0 = arma::log(1 - prob_prior);
   arma::vec lpr_prior_1 = arma::log(prob_prior);
   
   int tt = y.n_elem / k;
   
   bool tvp = false;
   if (a.n_elem / m > 1){
     tvp = true;
   }
   
   arma::sp_mat zz;
   if (tvp) {
     zz = sur_const_to_tvp(z, k, tt);
   }
   
   arma::mat a_mat;
   if (tvp) {
     a_mat = arma::reshape(a, m, tt); 
   } else {
     a_mat = a;
   }
   arma::mat AG = lambda * a_mat;
   arma::mat theta0 = AG;
   arma::mat theta1 = AG;
   arma::vec l0_res = arma::zeros<arma::vec>(k * tt);
   arma::vec l1_res = arma::zeros<arma::vec>(k * tt);
   //arma::vec g = arma::zeros<arma::vec>(1);
   
   double l0 = 0;
   double l1 = 0;
   double bayes = 0;
   double bayes_rand = 0;
   
   // If specified, limit the number of coefficients, where BVS is used
   arma::vec pos_res(m);
   for (arma::uword l = 0; l < m; l++) {
     pos_res(l) = l;
   }
   arma::uvec ex;
   if (include.isNotNull()) {
     ex = Rcpp::as<arma::uvec>(include) - 1;
     pos_res = pos_res.elem(ex);
   }
   pos_res = shuffle(pos_res);
   int pos_size = size(pos_res)(0);
   int var;
   
   // Start testing
   for (int j = 0; j < pos_size; j++){
     var = pos_res(j);
     //g = log(arma::randu<arma::vec>(1));
     //if (lambda(var, var) == 1 && g(0) >= lpr_prior_1(var)){continue;}
     //if (lambda(var, var) == 0 && g(0) >= lpr_prior_0(var)){continue;}
     //if ((lambda(var, var) == 1 && g(0) < lpr_prior_1(var)) || (lambda(var, var) == 0 && g(0) < lpr_prior_0(var))){
       theta0 = AG;
       theta1 = AG;
       theta1.row(var) = a_mat.row(var);
       if (tvp) {
         theta0.row(var) = arma::zeros<arma::mat>(1, tt);
         l0_res = y - zz * arma::vectorise(theta0);
         l1_res = y - zz * arma::vectorise(theta1);
       } else {
         theta0.row(var) = 0;
         l0_res = y - z * arma::vectorise(theta0);
         l1_res = y - z * arma::vectorise(theta1);
       }
       l0 = -.5 * arma::as_scalar(trans(l0_res) * sigma_i * l0_res) + arma::as_scalar(lpr_prior_0(var));
       l1 = -.5 * arma::as_scalar(trans(l1_res) * sigma_i * l1_res) + arma::as_scalar(lpr_prior_1(var));
       bayes = l1 - l0;
       bayes_rand = log(arma::as_scalar(arma::randu<arma::vec>(1)));
       if (bayes >= bayes_rand){
         lambda(var, var) = 1;
       } else {
         lambda(var, var) = 0;
       }
     //}
   }
   return lambda;
}

/*** R
# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model input, which uses BVS as variable selection algorithm
object <- create_bvarmodel(data = e1, p = 2, deterministic = "const",
                           varsel = "bvs")

# Add prior specifications, including the prior inclusion probabilities
object <- add_priors(object,
                     coef = list(v_i = 1, v_i_det = 1 / 10),
                     sigma = list(df = "k", scale = 1),
                     varsel = list(inprior = .1))

# Add initial values
object <- add_initial_values(object)

# Obtain data, initial values and priors
y <- matrix(t(object[["data"]][["train"]][["y"]]))
z <- object[["data"]][["train"]][["z"]] # Argument 'z' is taken dense
k <- object[["model"]][["k"]]
tt <- nrow(object[["data"]][["train"]][["y"]])
m <- ncol(z)
a <- object[["initial"]][["a"]]
prob_prior <- object[["priors"]][["a"]][["inprior"]]

# Arguments 'lambda' and 'sigma_i' have to be sparse
lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)

# Initial value of the inverse error covariance matrix
u <- matrix(y - z %*% a, k)
sigma_i <- Matrix::Matrix(kronecker(diag(1, tt), solve(tcrossprod(u) / tt)),
                          sparse = TRUE)

# Draw inclusion parameters
post_bvs(y, z, a, k, m, lambda, sigma_i, prob_prior)
*/