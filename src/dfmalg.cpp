// Important: this definition ensures Armadillo enables SuperLU. Use in R: Sys.setenv("PKG_LIBS"="-lsuperlu")
//#define ARMA_USE_SUPERLU 1

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.dfmalg)]]
Rcpp::List dfmalg(Rcpp::List object) {
  
  // Initialise variables
  Rcpp::List data = object["data"];
  
  arma::mat x = Rcpp::as<arma::mat>(data["X"]);
  arma::mat xvec = arma::vectorise(arma::trans(x));
  
  Rcpp::List model = object["model"];
  Rcpp::List factor = model["factors"];
  
  int tt = x.n_rows;
  int m = x.n_cols;
  int n = factor["number"];
  int p = factor["lags"];
  int n_lambda = (2 * m - n - 1) * n / 2;
  int n_phi =  n * n * p;
  
  // Priors ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  // Lambda
  Rcpp::List priors_lambda = priors["lambda"];
  arma::mat prior_lambda_Vi = Rcpp::as<arma::mat>(priors_lambda["v_i"]);
  // Sigma_u
  Rcpp::List priors_sigma_u = priors["sigma_u"];
  arma::vec prior_sigma_u_shape = priors_sigma_u["shape"];
  arma::vec prior_sigma_u_rate = priors_sigma_u["rate"];
  arma::vec post_sigma_u_shape = prior_sigma_u_shape + tt * 0.5;
  arma::vec post_sigma_u_scale = prior_sigma_u_rate * 0;
  
  // Phi
  Rcpp::List priors_phi;
  arma::mat prior_phi_mu;
  arma::mat prior_phi_Vi;
  if (p > 0) {
    priors_phi = priors["a"];
    prior_phi_mu = Rcpp::as<arma::mat>(priors_phi["mu"]);
    prior_phi_Vi = Rcpp::as<arma::mat>(priors_phi["v_i"]);
  }
  // Simga_v
  Rcpp::List priors_sigma_v = priors["sigma_v"];
  arma::vec prior_sigma_v_shape = priors_sigma_v["shape"];
  arma::vec prior_sigma_v_rate = priors_sigma_v["rate"];
  arma::vec post_sigma_v_shape = prior_sigma_v_shape + tt * 0.5;
  arma::vec post_sigma_v_scale = prior_sigma_v_rate * 0;
  
  // Initial values ----
  
  // Factor
  arma::sp_mat diag_tt = arma::speye<arma::sp_mat>(tt, tt);
  arma::sp_mat H_phi = arma::speye<arma::sp_mat>(tt * n, tt * n);
  arma::mat K_f;
  arma::mat K_f_dense = arma::mat(K_f);
  arma::mat f_hat, f, ff;
  
  // Lambda
  arma::mat K_l;
  int pos1, pos2, nlambda_i;
  arma::mat lambda_temp;
  arma::sp_mat lambda = arma::speye(m, n);
  arma::vec lambda_vec = arma::zeros<arma::vec>(n_lambda);
  arma::mat lambda_hat;
  
  // Sigma_u
  arma::mat u;
  arma::sp_mat sigma_i_u = arma::speye(m, m);
  sigma_i_u.diag() = 1 / arma::var(x);
  arma::vec sigma_u_temp(m);
  
  // Sigma_v
  arma::mat v;
  arma::sp_mat sigma_i_v = arma::speye(n, n);
  arma::vec sigma_v_temp(n);
  
  // Phi
  arma::vec fvec = arma::zeros<arma::vec>(n * tt);
  arma::mat K_phi, x_phi;
  arma::mat z_phi;
  arma::mat phi_hat;
  arma::mat phi;
  if (p > 0) {
    x_phi = arma::zeros<arma::mat>(tt, n * p);
    z_phi = arma::zeros<arma::mat>(tt * n, n * n * p);
    phi = arma::zeros<arma::mat>(n, n * p);
  }
  
  // Storage objects ----
  int iter = Rcpp::as<int>(model["iterations"]);
  int burnin = Rcpp::as<int>(model["burnin"]);
  int draws = iter + burnin;
  int pos_draw = 0;
  
  arma::mat draws_lambda = arma::zeros<arma::mat>(m * n, iter);
  arma::mat draws_phi;
  if (p > 0 ) {
    draws_phi = arma::zeros<arma::mat>(n_phi, iter); 
  }
  arma::mat draws_factor = arma::zeros<arma::mat>(n * tt, iter);
  arma::mat draws_sigma_u = arma::zeros<arma::mat>(m, iter);
  arma::mat draws_sigma_v = arma::mat(n, iter);
  
  // Gibbs sampler
  for (int draw = 0; draw < draws; draw++) {

    if (draw % 500 == 0) { // Check for user interuption ever now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Update H_phi
    if (p > 0) {
      phi = arma::reshape(phi, n, n * p);
      for (int i = 0; i < p; i++) {
        for (int j = 0; j < (tt - i - 1); j++) {
          H_phi.submat((j + 1) * n + i * n, j * n, (j + 1) * n + (i + 1) * n - 1, (j + 1) * n - 1) = -phi.cols(i * n, (i + 1) * n - 1);
        }
      }
    }
    
    // Draw factor
    K_f = arma::trans(H_phi) * arma::kron(diag_tt, sigma_i_v) * H_phi + arma::kron(diag_tt, arma::trans(lambda) * sigma_i_u * lambda);
    f_hat = arma::solve(K_f, arma::kron(diag_tt, arma::trans(lambda) * sigma_i_u) * xvec);
    K_f_dense = K_f;
    f = f_hat + arma::solve(chol(K_f_dense, "lower"), arma::randn(tt * n)); // MT x 1
    ff = arma::reshape(f, n, tt); // M x T
    fvec = arma::vectorise(ff);

    // Draw lambda
    int lambda_count = 0;
    for (int i = 1; i < m; i++) {
      pos1 = lambda_count;
      if (i <= (n - 1)) {
        nlambda_i = i;
        pos2 = lambda_count + i - 1;
        K_l = prior_lambda_Vi.submat(pos1, pos1, pos2, pos2) + ff.rows(0, i - 1) * ff.rows(0, i - 1).t() * sigma_i_u(i, i);
        lambda_hat = arma::solve(K_l, ff.rows(0, i - 1) * (x.col(i) - arma::trans(ff.row(i))) * sigma_i_u(i, i));
      } else {
        pos2 = lambda_count + n - 1;
        nlambda_i = n;
        K_l = prior_lambda_Vi.submat(pos1, pos1, pos2, pos2) + ff * ff.t() * sigma_i_u(i, i);
        lambda_hat = arma::solve(K_l, ff * x.col(i) * sigma_i_u(i, i));
      }
      lambda_temp = arma::trans(lambda_hat + arma::solve(chol(K_l, "lower"), arma::randn(nlambda_i)));
      lambda.submat(i, 0, i, nlambda_i - 1) = lambda_temp;
      lambda_count = pos2 + 1;
    }

    // Draw Sigma_u
    ff = arma::trans(ff);
    u = x - ff * arma::trans(lambda);
    post_sigma_u_scale = 1 / (prior_sigma_u_rate + arma::trans(arma::sum(arma::pow(u, 2))) * 0.5);
    for (int i = 0; i < m; i++) {
      sigma_i_u(i, i) = arma::randg<double>(arma::distr_param(post_sigma_u_shape(i), post_sigma_u_scale(i)));
    }
    
    // Draw Sigma_v
    v = arma::trans(arma::reshape(H_phi * fvec, n, tt));
    post_sigma_v_scale = 1 / (prior_sigma_v_rate + arma::trans(arma::sum(arma::pow(v, 2)) * 0.5));
    for (int i = 0; i < n; i++) {
      sigma_i_v(i, i) = arma::randg<double>(arma::distr_param(post_sigma_v_shape(i), post_sigma_v_scale(i)));
    }

    // Draw Phi
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        x_phi.submat(i + 1, n * i, tt - 1, n * (i + 1) - 1) = ff.rows(0, tt - 2 - i);
      }
      z_phi = arma::kron(x_phi, arma::eye(n, n));

      K_phi = prior_phi_Vi + arma::trans(z_phi) * arma::kron(diag_tt, sigma_i_v) * z_phi;
      phi_hat = arma::solve(K_phi, prior_phi_Vi * prior_phi_mu + arma::trans(z_phi) * arma::kron(diag_tt, sigma_i_v) * fvec);
      phi = phi_hat + arma::solve(chol(K_phi, "lower"), arma::randn(n * n * p));
    }

    // Store draws
    if (draw >= burnin) {
      pos_draw = draw - burnin;
      lambda_vec = arma::vectorise(lambda);
      draws_lambda.col(pos_draw) = lambda_vec;
      draws_factor.col(pos_draw) = fvec;
      sigma_u_temp = sigma_i_u.diag();
      draws_sigma_u.col(pos_draw) = 1 / sigma_u_temp;
      if (p > 0) {
        draws_phi.col(pos_draw) = phi;
      }
      sigma_v_temp = sigma_i_v.diag();
      draws_sigma_v.col(pos_draw) = 1 / sigma_v_temp;
    }
  } // End loop

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("lambda") = R_NilValue,
                                             Rcpp::Named("factor") = R_NilValue,
                                             Rcpp::Named("sigma_u") = R_NilValue,
                                             Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("sigma_v") = R_NilValue);

  posteriors["lambda"] = draws_lambda;
  posteriors["factor"] = draws_factor;
  // Changed notation after implementation
  posteriors["sigma_u"] = draws_sigma_u;
  if (p > 0) {
    posteriors["a"] = draws_phi;
  }
  posteriors["sigma_v"] = draws_sigma_v;

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posteriors") = posteriors);

  // Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") =  lambda_temp);
  // return result;
  
}