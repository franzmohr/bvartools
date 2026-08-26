// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export(.generate_forecasts)]]
arma::mat generate_forecasts(int k,
                        int p,
                        int h,
                        bool use_a,
                        arma::mat z,
                        arma::mat draws_a,
                        arma::mat draws_u_sigma_inv) {
  

  int draws = draws_u_sigma_inv.n_rows;
  
  arma::mat fcst = arma::mat(h * k, draws);
  const arma::mat diag_k = arma::eye<arma::mat>(k, k);
  arma::vec eigval;
  arma::mat eigvec;
  
  // Calculate forecasts
  for (int draw = 0; draw < draws; draw++)
  {
    for (int i = 0; i < h; i++)
    {
      if (use_a)
      {
        // Update z
        if (i > 0 && p > 0)
        {
          if (i < p)
          {
            z.submat(i * k, 0, (i + 1) * k - 1, i * k * k - 1) = arma::kron(arma::trans(fcst.submat(0, draw, i * k - 1, draw)), diag_k);
          }
          else
          {
            z.submat(i * k, 0, (i + 1) * k - 1, p * k * k - 1) = arma::kron(arma::trans(fcst.submat((i - p) * k, draw, i * k - 1, draw)), diag_k);
          }
        }
        // Update forecast
        fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = z.rows(i * k, (i + 1) * k - 1) * arma::trans(draws_a.row(draw));
      }
      
      // Add error
      arma::eig_sym(eigval, eigvec, arma::solve(arma::reshape(draws_u_sigma_inv.row(draw), k, k), diag_k));
      fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::trans(eigvec) * arma::randn(k);
    }
  }
  
  return arma::trans(fcst);
}