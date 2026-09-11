// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include "impact_matrix.h"

// [[Rcpp::export(.vardecomp)]]
arma::mat vardecomp(Rcpp::List A, int h, std::string type, int response) {
  
  // Get posterior draws and shock
  arma::mat a = Rcpp::as<arma::mat>(A["A"]);
  std::string oir ("oir");
  std::string sir ("sir");
  std::string gir ("gir");
  std::string sgir ("sgir");
  std::string custom ("custom");

  if (type != oir && type != sir && type != gir && type != sgir &&
      type != custom) {
    Rcpp::stop("vardecomp: unknown type \"%s\".", type);
  }

  // Collect information data
  int k = a.n_rows; // Number of endogenous variables
  const int lags = a.n_cols / k; // Lag order
  // Regressor columns the recursion can reach: lags beyond the horizon never
  // enter it. A count of columns, not a lag order, and zero on impact -- where
  // there is no recursion and nothing to slice.
  int n_use = a.n_cols;
  if (h < lags) {
    n_use = h * k;
  }
  
  arma::mat sigma = Rcpp::as<arma::mat>(A["Sigma"]);
  arma::mat a0i, sigma_mse, temp;
  arma::mat P = arma::eye<arma::mat>(k, k);
  
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"])));
    sigma_mse = sigma;
  }
  if (type == sir) {
    // The structural shocks carry the variances in Sigma, so the impulse matrix
    // is A_0^{-1} times a factor of Sigma and the forecast error covariance is
    // A_0^{-1} Sigma A_0^{-1}'. Dropping Sigma here would decompose the variance
    // as if every structural shock had unit variance, which shifts weight to
    // whichever shock has the largest loading in A_0.
    a0i = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k));
    P = a0i * arma::trans(arma::chol(sigma));
    sigma_mse = a0i * sigma * arma::trans(a0i);
  }
  if (type == gir) {
    P = sigma;
    sigma_mse = sigma;
  }
  if (type == sgir) {
    a0i = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k));
    P = a0i * sigma;
    sigma_mse = a0i * sigma * arma::trans(a0i);
  }
  if (type == custom) {
    // A caller supplied impact matrix relabels the shocks but says nothing
    // about the model, so the forecast error covariance stays Sigma. The
    // shares then add up across shocks exactly when P P' = Sigma -- true of a
    // rotation of the Cholesky factor, and not of an arbitrary matrix. Nothing
    // here enforces it: a decomposition that does not sum to one is the honest
    // report of an impact matrix that does not factorise Sigma.
    P = bvartools_impact_matrix(A, k, "fevd");
    sigma_mse = sigma;
  }
  
  double sigmajj = 1;
  if (type == gir || type == sgir) {
    sigmajj = sqrt(arma::as_scalar(sigma(response - 1, response - 1)));
  }
  
  // Matrix of coefficients. Only needed once the recursion below runs: on
  // impact the decomposition rests on Phi_0 = I alone.
  arma::mat A_temp;
  if (h > 0) {
    A_temp = arma::zeros<arma::mat>(k, h * k);
    if (n_use > 0) {
      A_temp.cols(0, n_use - 1) = a.cols(0, n_use - 1);
    }
  }
  
  // Generate output object
  arma::mat result = arma::zeros<arma::mat>(h + 1, k);
  arma::mat numerator = result;
  arma::vec mse = arma::zeros<arma::vec>(h + 1);
  arma::mat ejt = arma::zeros<arma::mat>(1, k);
  ejt(0, response - 1) = 1;
  
  // Generate FEIRF
  arma::mat phi = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi_temp = arma::eye<arma::mat>(k, k);
  phi.rows(0, k - 1) = phi_temp; // Set first element of FEIRF to identity matrix
  
  // Time = 0
  numerator.row(0) = arma::square(ejt * phi_temp * P);
  mse(0) = arma::as_scalar(ejt * phi_temp * sigma_mse * arma::trans(phi_temp) * arma::trans(ejt));
  result.row(0) =  numerator.row(0) / sigmajj / mse(0);
  
  for (int i = 1; i <= h; i++) {
    // FEIR
    phi_temp.zeros(); // Reset phi_temp
    for (int j = 1; j <= i; j++) {
      phi_temp = phi_temp + phi.rows((i - j) * k, (i - j + 1) * k - 1) * A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi.rows(i * k, (i + 1) * k - 1) = phi_temp ; // Update FEIRF
    
    // Generate GFEVD
    numerator.row(i) = numerator.row(i - 1) + arma::square(ejt * phi_temp * P);
    mse(i) = mse(i - 1) +  arma::as_scalar(ejt * phi_temp * sigma_mse * arma::trans(phi_temp) * arma::trans(ejt));
    result.row(i) =  numerator.row(i) / sigmajj / mse(i); 
  }
  
  return result;
}