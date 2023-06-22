// Important: this definition ensures Armadillo enables SuperLU. Use in R: Sys.setenv("PKG_LIBS"="-lsuperlu")
// #define ARMA_USE_SUPERLU 1

#include "../inst/include/bvartools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.bvartvpalg)]]
Rcpp::List bvartvpalg(Rcpp::List object) {
  
  // Model data and helper object
  Rcpp::List data = object["data"];
  arma::mat y = arma::trans(Rcpp::as<arma::mat>(data["Y"]));
  const arma::mat yvec = arma::vectorise(y);
  Rcpp::Nullable<Rcpp::List> z_test = data["SUR"];
  arma::mat z;
  int n_tot = 0;
  if (z_test.isNotNull()) {
    z = Rcpp::as<arma::mat>(data["SUR"]);
    n_tot = z.n_cols;
  }
  const bool use_a = n_tot > 0;
  
  // Model information
  Rcpp::List model = object["model"];
  Rcpp::CharacterVector model_names = model.names();
  Rcpp::List endogen = model["endogen"];
  Rcpp::CharacterVector endo_names = Rcpp::as<Rcpp::CharacterVector>(endogen["variables"]);
  
  const int tt = y.n_cols; // Number of periods
  const int k = y.n_rows; // Number of endogeous variable in the measurement equation
  const int p = Rcpp::as<int>(endogen["lags"]); // Lags of endogenous coefficiens
  const int n_a = k * k * p;
  int m = 0;
  int s = 0;
  int n_b = 0;
  int n = 0;
  int n_c = 0;
  const int n_sigma = k * k;
  const arma::mat diag_k = arma::eye<arma::mat>(k, k); // K diag matrix
  const arma::sp_mat diag_tt = arma::speye<arma::sp_mat>(tt, tt); // T diag matrix
  const arma::vec vec_tt = arma::ones<arma::vec>(tt); // T vector
  const bool sv = Rcpp::as<bool>(model["sv"]);
  const bool structural = Rcpp::as<bool>(model["structural"]);
  int n_a0 = 0;
  int n_psi = 0;
  if (structural) {
    n_a0 = k * (k - 1) / 2;
  }
  
  bool covar = false;
  bool bvs = false;
  bool psi_bvs = false;
  
  // Exogenous variables
  Rcpp::CharacterVector exogen_names;
  Rcpp::List exogen;
  if (std::find(model_names.begin(), model_names.end(), "exogen") != model_names.end()) {
    exogen = model["exogen"];
    exogen_names = Rcpp::as<Rcpp::CharacterVector>(exogen["variables"]);
    m = exogen_names.length();
    s = Rcpp::as<int>(exogen["lags"]);
    n_b = k * m * (s + 1);
  }
  
  // Deterministic terms
  Rcpp::CharacterVector det_names;
  if (std::find(model_names.begin(), model_names.end(), "deterministic") != model_names.end()) {
    det_names = Rcpp::as<Rcpp::CharacterVector>(model["deterministic"]);
    n = det_names.length();
    n_c = k * n;
  }
  
  Rcpp::List initial = object["initial"];
  Rcpp::List init_coeffs;
  Rcpp::List init_psi;
  Rcpp::List init_sigma = initial["sigma"];
  
  // Priors & initial values ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  
  // Priors - Coefficients
  Rcpp::List priors_coefficients, a_prior_varsel;
  Rcpp::CharacterVector prcoeff_names;
  arma::vec a, a_init, a_init_post_mu, a_lag, a_post_mu, a_sigma_v_post_shape, a_sigma_v_prior_rate, a_sigma_v_post_scale;
  arma::vec a_bvs_lprior_0, a_bvs_lprior_1, a_bvs_l0_res, a_bvs_l1_res, a_randu, a_prior_incl, a_varsel_include, a_varsel_include_draw;
  arma::mat a_AG, a_init_prior_mu, a_init_prior_vi, a_init_post_v, a_mat, a_post_v, a_theta0, a_theta1, a_v;
  arma::sp_mat a_hh, a_lambda, a_sigma_v_i, zz, a_zzss_i, a_h_sigmav_i_h, zz_bvs;
  int a_varsel_n, a_varsel_pos;
  double a_l0, a_l1, a_bayes, a_bayes_rand;
  if (use_a) {
    priors_coefficients = priors["coefficients"];
    prcoeff_names = priors_coefficients.names();
    zz = arma::zeros<arma::sp_mat>(tt * k, tt * n_tot); // Final data matrix
    for (int i = 0; i < tt; i++) {
      zz.submat(i * k, i * n_tot, (i + 1) * k - 1, (i + 1) * n_tot - 1) = z.rows(i * k, (i + 1) * k - 1);
    } 
    // Measurement
    a_init_prior_mu = Rcpp::as<arma::mat>(priors_coefficients["mu"]);
    a_init_prior_vi = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
    a_init_post_mu = a_init_prior_mu * 0;
    // State
    a_sigma_v_post_shape = Rcpp::as<arma::vec>(priors_coefficients["shape"]) + 0.5 * tt;
    a_sigma_v_prior_rate = Rcpp::as<arma::vec>(priors_coefficients["rate"]);
    
    a_post_mu = arma::zeros<arma::vec>(n_tot * tt);
    a_post_v = arma::zeros<arma::mat>(tt * n_tot, tt * n_tot);
    a = arma::zeros<arma::vec>(n_tot * tt);
    a_lag = a;
    a_sigma_v_i = arma::speye<arma::sp_mat>(n_tot, n_tot);
    a_sigma_v_i.diag() = 1 / a_sigma_v_prior_rate;
    
    init_coeffs = initial["coefficients"];
    a_init = Rcpp::as<arma::vec>(init_coeffs["draw"]);
    a_hh = arma::speye<arma::sp_mat>(n_tot * tt, n_tot * tt); // TVP H matrix
    a_hh.diag(-n_tot) = -arma::ones<arma::vec>(n_tot * (tt - 1));

    // Priors - Coefficients - BVS    
    if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "bvs") != prcoeff_names.end()) {
      bvs = true;
      a_prior_varsel = priors_coefficients["bvs"];
      a_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
      a_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
      a_varsel_include = Rcpp::as<arma::vec>(a_prior_varsel["include"]) - 1;
      a_varsel_n = size(a_varsel_include)(0);
      a_lambda = arma::eye<arma::sp_mat>(n_tot, n_tot);
      a_l0 = 0;
      a_l1 = 0;
      a_bayes = 0;
      a_bayes_rand = 0;
      zz_bvs = zz;
    }
  }
  
  // Priors - Covar coefficients
  Rcpp::List psi_priors, psi_prior_varsel;
  Rcpp::CharacterVector psi_priors_names;
  arma::vec psi_prior_incl, psi_tau0, psi_tau1, psi_tau0sq, psi_tau1sq, psi_bvs_lprior_0, psi_bvs_lprior_1;
  arma::vec psi_sigma_v_post_scale, psi_sigma_v_post_shape, psi_sigma_v_prior_rate;
  arma::mat psi_init_prior_mu, psi_init_prior_vi;
  
  if (std::find(priors_names.begin(), priors_names.end(), "psi") != priors_names.end()) {
    covar = true;
    psi_priors = priors["psi"];
    psi_priors_names = psi_priors.names();
    psi_init_prior_mu = Rcpp::as<arma::mat>(psi_priors["mu"]);
    psi_init_prior_vi = Rcpp::as<arma::mat>(psi_priors["v_i"]);
    psi_sigma_v_post_shape = Rcpp::as<arma::vec>(psi_priors["shape"]) + 0.5 * tt;
    psi_sigma_v_prior_rate = Rcpp::as<arma::vec>(psi_priors["rate"]);
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "bvs") != psi_priors_names.end()) {
      psi_bvs = true;
      psi_prior_varsel = psi_priors["bvs"];
      psi_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
      psi_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
    }
  }
  // Initial values
  int psi_varsel_n, psi_varsel_pos;
  double psi_bayes, psi_bayes_rand, psi_l0, psi_l1;
  arma::vec psi, psi_init, psi_init_post_mu, psi_lag, psi_post_incl, psi_post_mu, psi_randu, psi_theta0_res, psi_theta1_res, psi_varsel_include, psi_varsel_include_draw, psi_u0, psi_u1, psi_y;
  arma::mat diag_Psi, psi_AG, psi_init_post_v, psi_mat, psi_post_v, psi_theta0, psi_theta1, psi_v, psi_z_bvs;
  arma::sp_mat diag_covar_omega_i, Psi, psi_hh, psi_h_sigmav_i_h, psi_lambda, psi_sigma_v_i, psi_z, psi_zzss_i;
  if (covar) {
    n_psi = k * (k - 1) / 2;
    Psi = arma::speye<arma::sp_mat>(k * tt, k * tt);
    Rcpp::List init_psi = initial["psi"];
    psi_init = Rcpp::as<arma::vec>(init_psi["draw"]);
    psi_lag = arma::zeros<arma::vec>(n_psi * tt);
    psi_z = arma::zeros<arma::mat>((k - 1) * tt, n_psi * tt);
    psi_hh = arma::speye<arma::sp_mat>(n_psi * tt, n_psi * tt);
    psi_hh.diag(-n_psi) = -arma::ones<arma::vec>(n_psi * (tt - 1));
    diag_covar_omega_i = arma::zeros<arma::sp_mat>(tt * (k - 1), tt * (k - 1));
    psi_sigma_v_i = arma::speye<arma::sp_mat>(n_psi, n_psi);
    psi_sigma_v_i.diag() = 1 / psi_sigma_v_prior_rate;
    if (psi_bvs) {
      psi_varsel_include = Rcpp::as<arma::vec>(psi_prior_varsel["include"]) - 1;
      psi_varsel_n = size(psi_varsel_include)(0);
      psi_lambda = arma::speye<arma::sp_mat>(n_psi, n_psi);
      psi_l0 = 0;
      psi_l1 = 0;
      psi_bayes = 0;
      psi_bayes_rand = 0;
    }
  }
  
  // Priors - Measurement errors
  Rcpp::List sigma_pr = priors["sigma"];
  Rcpp::CharacterVector sigma_names = sigma_pr.names();
  double sigma_post_df;
  arma::vec sigma_post_shape, sigma_prior_rate, sigma_prior_mu;
  arma::mat sigma_prior_scale, sigma_prior_vi;
  bool use_gamma = false;
  if (sv) {
    sigma_prior_mu = Rcpp::as<arma::vec>(sigma_pr["mu"]);
    sigma_prior_vi = Rcpp::as<arma::mat>(sigma_pr["v_i"]);
    sigma_post_shape = Rcpp::as<arma::vec>(sigma_pr["shape"]) + 0.5 * tt;
    sigma_prior_rate = Rcpp::as<arma::vec>(sigma_pr["rate"]);
  } else {
    if (std::find(sigma_names.begin(), sigma_names.end(), "df") != sigma_names.end()) {
      sigma_post_df = Rcpp::as<double>(sigma_pr["df"]) + tt;
      sigma_prior_scale = Rcpp::as<arma::mat>(sigma_pr["scale"]);
    }
    if (std::find(sigma_names.begin(), sigma_names.end(), "shape") != sigma_names.end()) {
      use_gamma = true;
      sigma_post_shape = Rcpp::as<arma::vec>(sigma_pr["shape"]) + 0.5 * tt;
      sigma_prior_rate = Rcpp::as<arma::vec>(sigma_pr["rate"]);
    }
  }
  // Initial values
  arma::vec h_init, sigma_h, u_vec, sigma_post_scale;
  arma::mat h_init_post_v, sigma_h_i, diag_sigma_i_temp;
  arma::vec h_init_post_mu;
  arma::mat h, h_lag, sse;
  arma::mat u = y * 0;
  arma::sp_mat diag_omega_i, sigma_u_i, omega_i;
  arma::sp_mat diag_sigma_u_i = arma::zeros<arma::sp_mat>(k * tt, k * tt);
  if (sv) {
    h = Rcpp::as<arma::mat>(init_sigma["h"]);
    h_lag = h * 0;
    sigma_h = Rcpp::as<arma::vec>(init_sigma["sigma_h"]);
    h_init = arma::vectorise(h.row(0));
    sigma_u_i = arma::diagmat(1 / exp(h_init));
  } else {
    omega_i = Rcpp::as<arma::mat>(init_sigma["sigma_i"]);
    sigma_u_i = omega_i;
  }
  diag_sigma_u_i.diag() = arma::repmat(sigma_u_i.diag(), tt, 1);
  if (covar || sv) {
    diag_omega_i = diag_sigma_u_i;
  }
  
  // Storage objects
  int iter = Rcpp::as<int>(model["iterations"]);
  int burnin = Rcpp::as<int>(model["burnin"]);
  int draws = iter + burnin;
  int pos_draw;
  const int a_pos_start = 0;
  const int a_pos_end = n_a - 1;
  const int b_pos_start = n_a;
  const int b_pos_end = n_a + n_b - 1;
  const int c_pos_start = n_a + n_b;
  const int c_pos_end = n_a + n_b + n_c - 1;
  const int a0_pos_start = n_a + n_b + n_c;
  const int a0_pos_end = n_a + n_b + n_c + n_a0 - 1;
  
  arma::mat draws_a0 = arma::zeros<arma::mat>(n_a0 * tt, iter);
  arma::mat draws_sigma_a0 = arma::zeros<arma::mat>(n_a0, iter);
  arma::mat draws_a = arma::zeros<arma::mat>(n_a * tt, iter);
  arma::mat draws_sigma_a = arma::zeros<arma::mat>(n_a, iter);
  arma::mat draws_b = arma::zeros<arma::mat>(n_b * tt, iter);
  arma::mat draws_sigma_b = arma::zeros<arma::mat>(n_b, iter);
  arma::mat draws_c = arma::zeros<arma::mat>(n_c * tt, iter);
  arma::mat draws_sigma_c = arma::zeros<arma::mat>(n_c, iter);
  arma::mat draws_sigma_u, draws_sigma_sigma;
  if (sv || covar) {
    draws_sigma_u = arma::zeros<arma::mat>(k * k * tt, iter);
  } else {
    draws_sigma_u = arma::zeros<arma::mat>(k * k, iter);
  }
  if (sv) {
    draws_sigma_sigma = arma::zeros<arma::mat>(k * k, iter);
  }
  
  arma::vec a_lambda_vec, psi_lambda_vec;
  arma::mat draws_lambda_a0, draws_lambda_a, draws_lambda_b, draws_lambda_c;
  if (bvs) {
    if (structural) {
      draws_lambda_a0 = arma::zeros<arma::mat>(n_a0, iter); 
    }
    draws_lambda_a = arma::zeros<arma::mat>(n_a, iter);
    draws_lambda_b = arma::zeros<arma::mat>(n_b, iter);
    draws_lambda_c = arma::zeros<arma::mat>(n_c, iter);
  }
  if (covar && psi_bvs) {
    draws_lambda_a0 = arma::zeros<arma::mat>(n_psi, iter);  
  }
  
  // Start Gibbs sampler
  for (int draw = 0; draw < draws; draw++) {

    if (draw % 20 == 0) { // Check for user interuption ever now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Draw a
    if (use_a) {
      if (bvs) {
        zz = zz_bvs * arma::kron(diag_tt, a_lambda);
      }
      a_zzss_i = arma::trans(zz) * diag_sigma_u_i;
      a_h_sigmav_i_h = arma::trans(a_hh) * arma::kron(diag_tt, a_sigma_v_i) * a_hh;
      a_post_v = a_h_sigmav_i_h + a_zzss_i * zz;
      a_post_mu = arma::solve(a_post_v, a_h_sigmav_i_h * arma::kron(vec_tt, a_init) + a_zzss_i * yvec);
      a = a_post_mu + arma::solve(arma::chol(a_post_v), arma::randn<arma::vec>(n_tot * tt));
      
      // BVS
      if (bvs) {
        zz = zz_bvs;
        a_mat = arma::reshape(a, n_tot, tt);
        a_AG = a_lambda * a_mat;
        a_varsel_include_draw = shuffle(a_varsel_include); // Reorder positions of variable selection
        for (int j = 0; j < a_varsel_n; j++){
          a_varsel_pos = a_varsel_include_draw(j);
          a_randu = arma::log(arma::randu<arma::vec>(1));
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) >= a_bvs_lprior_1(a_varsel_pos)){continue;}
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) >= a_bvs_lprior_0(a_varsel_pos)){continue;}
          if ((a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) < a_bvs_lprior_1(a_varsel_pos)) || (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) < a_bvs_lprior_0(a_varsel_pos))){
            a_theta0 = a_AG;
            a_theta1 = a_AG;
            a_theta0.row(a_varsel_pos) = arma::zeros<arma::mat>(1, tt);
            a_theta1.row(a_varsel_pos) = a_mat.row(a_varsel_pos);
            a_bvs_l0_res = yvec - zz * arma::vectorise(a_theta0);
            a_bvs_l1_res = yvec - zz * arma::vectorise(a_theta1);
            a_l0 = -arma::as_scalar(trans(a_bvs_l0_res) * diag_sigma_u_i * a_bvs_l0_res) * 0.5 + arma::as_scalar(a_bvs_lprior_0(a_varsel_pos));
            a_l1 = -arma::as_scalar(trans(a_bvs_l1_res) * diag_sigma_u_i * a_bvs_l1_res) * 0.5 + arma::as_scalar(a_bvs_lprior_1(a_varsel_pos));
            a_bayes = a_l1 - a_l0;
            a_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (a_bayes >= a_bayes_rand){
              a_lambda(a_varsel_pos, a_varsel_pos) = 1;
            } else {
              a_lambda(a_varsel_pos, a_varsel_pos) = 0;
            }
          }
        }
        a = arma::vectorise(a_lambda * a_mat);
        a_lambda_vec = a_lambda.diag();
      }
      
      u_vec = yvec - zz * a;
    } else {
      u_vec = yvec;
    }
    u = arma::reshape(u_vec, k, tt);
    
    // Covariances
    if (covar) {
      
      // Prepare data
      psi_y = arma::vectorise(u.rows(1, k - 1));
      for (int i = 1; i < k; i++) {
        for (int j = 0; j < tt; j++) {
          psi_z.submat(j * (k - 1) + i - 1,
                       j * n_psi + i * (i - 1) / 2,
                       j * (k - 1) + i - 1,
                       j * n_psi + (i + 1) * i / 2 - 1) = -arma::trans(u.submat(0, j, i - 1, j));
          
          diag_covar_omega_i(j * (k - 1) + i - 1, j * (k - 1) + i - 1) = diag_omega_i(j * k + i, j * k + i);
        }
      }
      
      if (psi_bvs) {
        psi_z_bvs = psi_z;
        psi_z = psi_z_bvs * arma::kron(diag_tt, psi_lambda);
      }
      
      psi_zzss_i = arma::trans(psi_z) * diag_covar_omega_i;
      psi_h_sigmav_i_h = arma::trans(psi_hh) * arma::kron(diag_tt, psi_sigma_v_i) * psi_hh;
      psi_post_v = psi_h_sigmav_i_h + psi_zzss_i * psi_z;
      psi_post_mu = arma::solve(psi_post_v, psi_h_sigmav_i_h * arma::kron(vec_tt, psi_init) + psi_zzss_i * psi_y);
      psi = psi_post_mu + arma::solve(arma::chol(psi_post_v), arma::randn(n_psi * tt));
      
      if (psi_bvs) {
        
        // Reorder positions of variable selection
        psi_varsel_include_draw = shuffle(psi_varsel_include);
        
        psi_z = psi_z_bvs;
        psi_mat = arma::reshape(psi, n_psi, tt);
        psi_AG = psi_lambda * psi_mat;
        for (int j = 0; j < psi_varsel_n; j++){
          psi_varsel_pos = psi_varsel_include_draw(j);
          psi_randu = arma::log(arma::randu<arma::vec>(1));
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) >= psi_bvs_lprior_1(psi_varsel_pos)){continue;}
          if (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) >= psi_bvs_lprior_0(psi_varsel_pos)){continue;}
          if ((psi_lambda(psi_varsel_pos, psi_varsel_pos) == 1 && psi_randu(0) < psi_bvs_lprior_1(psi_varsel_pos)) || (psi_lambda(psi_varsel_pos, psi_varsel_pos) == 0 && psi_randu(0) < psi_bvs_lprior_0(psi_varsel_pos))){
            psi_theta0 = psi_AG;
            psi_theta1 = psi_AG;
            psi_theta0.row(psi_varsel_pos) = 0;
            psi_theta1.row(psi_varsel_pos) = psi.row(psi_varsel_pos);
            psi_theta0_res = psi_y - psi_z * psi_theta0;
            psi_theta1_res = psi_y - psi_z * psi_theta1;
            psi_l0 = -arma::as_scalar(trans(psi_theta0_res) * diag_covar_omega_i * psi_theta0_res) / 2 + arma::as_scalar(psi_bvs_lprior_0(psi_varsel_pos));
            psi_l1 = -arma::as_scalar(trans(psi_theta1_res) * diag_covar_omega_i * psi_theta1_res) / 2 + arma::as_scalar(psi_bvs_lprior_1(psi_varsel_pos));
            psi_bayes = psi_l1 - psi_l0;
            psi_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (psi_bayes >= psi_bayes_rand){
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 1;
            } else {
              psi_lambda(psi_varsel_pos, psi_varsel_pos) = 0;
            }
          }
        }
        psi = arma::vectorise(psi_lambda * psi_mat);
        psi_lambda_vec = psi_lambda.diag();
      }

      for (int j = 0; j < tt; j ++) {
        for (int i = 1; i < k; i++) {
          Psi.submat((k * j) + i, k * j, (k * j) + i, (k * j) + i - 1) = arma::trans(psi.subvec((n_psi * j) + i * (i - 1) / 2, (n_psi * j) + (i + 1) * i / 2 - 1));
        }
      }
      
      u = arma::reshape(Psi * u_vec, k, tt);
    }
    
    if (sv) {

      // Draw variances
      for (int i = 0; i < k; i++) {
        h.col(i) = bvartools::stoch_vol(u.row(i).t(), h.col(i), sigma_h(i), h_init(i));
      }
      diag_omega_i.diag() = 1 / exp(arma::vectorise(h.t()));
      if (covar) {
        diag_sigma_u_i = arma::trans(Psi) * diag_omega_i * Psi;
      } else {
        diag_sigma_u_i = diag_omega_i;
      }

      // Draw sigma_h
      h_lag.row(0) = h_init.t();
      h_lag.rows(1, tt - 1) = h.rows(0, tt - 2);
      h_lag = h - h_lag;
      sigma_post_scale = 1 / (sigma_prior_rate + arma::trans(arma::sum(arma::pow(h_lag, 2))) * 0.5);
      for (int i = 0; i < k; i++) {
        sigma_h(i) = 1 / arma::randg<double>(arma::distr_param(sigma_post_shape(i), sigma_post_scale(i)));
      }

      // Draw h_init
      sigma_h_i = arma::diagmat(1 / sigma_h);
      h_init_post_v = sigma_prior_vi + sigma_h_i;
      h_init_post_mu = arma::solve(h_init_post_v, sigma_prior_vi * sigma_prior_mu + sigma_h_i * h.row(0).t());
      h_init = h_init_post_mu + arma::solve(arma::chol(h_init_post_v), arma::randn(k));

    } else {

      if (use_gamma) {

        sse = u * u.t();
        for (int i = 0; i < k; i++) {
          omega_i(i, i) = arma::randg<double>(arma::distr_param(sigma_post_shape(i), 1 / arma::as_scalar(sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }
        diag_omega_i = arma::kron(diag_tt, omega_i);
        if (covar) {
          diag_sigma_u_i = arma::trans(Psi) * diag_omega_i * Psi;
        } else {
          sigma_u_i = omega_i;
          diag_sigma_u_i = diag_omega_i;
        }

      } else {
        sigma_u_i = arma::wishrnd(arma::solve(sigma_prior_scale + u * u.t(), diag_k), sigma_post_df);
        diag_sigma_u_i = arma::kron(diag_tt, sigma_u_i);
      }
    }

    if (use_a) {
      // Draw sigma_v_i
      a_lag.subvec(0, n_tot - 1) = a_init;
      a_lag.subvec(n_tot, tt * n_tot - 1) = a.subvec(0, (tt - 1) * n_tot - 1);
      a_v = arma::trans(arma::reshape(a - a_lag, n_tot, tt));
      a_sigma_v_post_scale = 1 / (a_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(a_v, 2))) * 0.5);
      for (int i = 0; i < n_tot; i++) {
        a_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(a_sigma_v_post_shape(i), a_sigma_v_post_scale(i)));
      }
      
      // Draw initial state of a
      a_init_post_v = a_init_prior_vi + a_sigma_v_i;
      a_init_post_mu = arma::solve(a_init_post_v, a_init_prior_vi * a_init_prior_mu + a_sigma_v_i * a.subvec(0, n_tot - 1));
      a_init = a_init_post_mu + arma::solve(arma::chol(a_init_post_v), arma::randn(n_tot)); 
    }
    
    if (covar) {
      // Draw sigma_v_i
      psi_lag.subvec(0, n_psi - 1) = psi_init;
      psi_lag.subvec(n_psi, tt * n_psi - 1) = psi.subvec(0, (tt - 1) * n_psi - 1);
      psi_v = arma::trans(arma::reshape(psi - psi_lag, n_psi, tt));
      psi_sigma_v_post_scale = 1 / (psi_sigma_v_prior_rate + arma::vectorise(arma::sum(arma::pow(psi_v, 2))) * 0.5);
      for (int i = 0; i < n_psi; i++) {
        psi_sigma_v_i(i, i) = arma::randg<double>(arma::distr_param(psi_sigma_v_post_shape(i), psi_sigma_v_post_scale(i)));
      }

      // Draw initial state of a
      psi_init_post_v = psi_init_prior_vi + psi_sigma_v_i;
      psi_init_post_mu = arma::solve(psi_init_post_v, psi_init_prior_vi * psi_init_prior_mu + psi_sigma_v_i * psi.subvec(0, n_psi - 1));
      psi_init = psi_init_post_mu + arma::solve(arma::chol(psi_init_post_v), arma::randn(n_psi));
    }

    // Store draws
    if (draw >= burnin) {

      pos_draw = draw - burnin;

      if (sv || covar) {
        for (int i = 0; i < tt; i ++) {
          draws_sigma_u.submat(i * n_sigma, pos_draw, (i + 1) * n_sigma - 1, pos_draw) = arma::vectorise(arma::solve(arma::mat(diag_sigma_u_i.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1)), diag_k));
        }
        if (sv) {
          draws_sigma_sigma.col(pos_draw) = arma::vectorise(arma::diagmat(sigma_h));
        }
      } else {
        draws_sigma_u.col(pos_draw) = arma::vectorise(arma::solve(arma::mat(sigma_u_i), diag_k)); 
      }

      if (psi_bvs) {
        draws_lambda_a0.col(pos_draw) = psi_lambda_vec;
      }

      a_mat = arma::reshape(a, n_tot, tt);

      if (n_a > 0) {
        draws_a.col(pos_draw) = arma::vectorise(a_mat.rows(a_pos_start, a_pos_end));
        draws_sigma_a.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(a_pos_start, a_pos_start, a_pos_end, a_pos_end)).diag();
        if (bvs) {
          draws_lambda_a.col(pos_draw) = a_lambda_vec.subvec(a_pos_start, a_pos_end);
        }
      }
      if (n_b > 0) {
        draws_b.col(pos_draw) = arma::vectorise(a_mat.rows(b_pos_start, b_pos_end));
        draws_sigma_b.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(b_pos_start, b_pos_start, b_pos_end, b_pos_end)).diag();
        if (bvs) {
          draws_lambda_b.col(pos_draw) = a_lambda_vec.subvec(b_pos_start, b_pos_end);
        }
      }
      if (n_c > 0) {
        draws_c.col(pos_draw) = arma::vectorise(a_mat.rows(c_pos_start, c_pos_end));
        draws_sigma_c.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(c_pos_start, c_pos_start, c_pos_end, c_pos_end)).diag();
        if (bvs) {
          draws_lambda_c.col(pos_draw) = a_lambda_vec.subvec(c_pos_start, c_pos_end);
        }
      }
      if (structural) {
        draws_a0.col(pos_draw) = arma::vectorise(a_mat.rows(a0_pos_start, a0_pos_end));
        draws_sigma_a0.col(pos_draw) = 1 / arma::mat(a_sigma_v_i.submat(a0_pos_start, a0_pos_start, a0_pos_end, a0_pos_end)).diag();
        if (bvs) {
          draws_lambda_a0.col(pos_draw) = a_lambda_vec.subvec(a0_pos_start, a0_pos_end);
        }
      }
    } // End storage condition

  } // End loop

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("b") = R_NilValue,
                                             Rcpp::Named("c") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);

  if (n_a > 0) {
    if (bvs) {
      posteriors["a"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a,
                                                     Rcpp::Named("sigma") = draws_sigma_a,
                                                     Rcpp::Named("lambda") = draws_lambda_a));
    } else {
      posteriors["a"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a,
                                                     Rcpp::Named("sigma") = draws_sigma_a));
    }
  }

  if (n_b > 0) {
    if (bvs) {
      posteriors["b"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_b,
                                                     Rcpp::Named("sigma") = draws_sigma_b,
                                                     Rcpp::Named("lambda") = draws_lambda_b));
    } else {
      posteriors["b"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_b,
                                                     Rcpp::Named("sigma") = draws_sigma_b));
    }
  }

  if (n_c > 0) {
    if (bvs) {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c,
                                                     Rcpp::Named("sigma") = draws_sigma_c,
                                                     Rcpp::Named("lambda") = draws_lambda_c));
    } else {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c,
                                                     Rcpp::Named("sigma") = draws_sigma_c));
    }
  }

  if (structural) {
    if (bvs) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                       Rcpp::Named("sigma") = draws_sigma_a0,
                                                       Rcpp::Named("lambda") = draws_lambda_a0));
    } else {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                       Rcpp::Named("sigma") = draws_sigma_a0));
    }
  }

  if (psi_bvs) {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    }
  } else {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma_u));
    }
  }

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posteriors") = posteriors);
  
  // return Rcpp::List::create(Rcpp::Named("test") = a);
  
}


