#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.vardecomp)]]
arma::mat vardecomp(Rcpp::List A, int h, std::string type, int response) {
  
  // Get posterior draws and shock
  arma::mat a = Rcpp::as<arma::mat>(A["A"]);
  
  std::string oir ("oir");
  std::string sir ("sir");
  std::string gir ("gir");
  std::string sgir ("sgir");
  
  // Collect information data
  int k = a.n_rows; // Number of endogenous variables
  int p = a.n_cols / k; // Lag order
  // If horizon is smaller than lag order...
  if (h < p) {
    p = h * k; // ...use reduced lag order for IRF
  } else {
    p = a.n_cols; // ...if not, use total number of regressors
  }
  
  arma::mat sigma = Rcpp::as<arma::mat>(A["Sigma"]);
  arma::mat a0i, sigma_mse, temp;
  arma::mat P = arma::eye<arma::mat>(k, k);
  
  if (type == oir) {
    P = arma::trans(arma::chol(Rcpp::as<arma::mat>(A["Sigma"])));
    sigma_mse = sigma;
  }
  if (type == sir) {
    a0i = arma::solve(Rcpp::as<arma::mat>(A["A0"]), arma::eye<arma::mat>(k, k));
    P = a0i;
    sigma_mse = a0i * arma::trans(a0i);
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
  
  double sigmajj = 1;
  if (type == gir || type == sgir) {
    sigmajj = sqrt(arma::as_scalar(sigma(response - 1, response - 1)));
  }
  
  // Matrix of coefficients
  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, p - 1) = a.cols(0, p - 1); 
  
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