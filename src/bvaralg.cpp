#include "../inst/include/bvartools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.bvaralg)]]
Rcpp::List bvaralg(Rcpp::List object) {
  
  // Initialise variables
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
  
  // Define useful variables
  const int tt = y.n_cols;
  const int k = y.n_rows;
  const int p = Rcpp::as<int>(endogen["lags"]);
  const int n_a = k * k * p;
  int m = 0;
  int s = 0;
  int n_b = 0;
  int n = 0;
  int n_c = 0;
  int n_psi = 0;
  const int n_sigma = k * k;
  const arma::mat diag_k = arma::eye<arma::mat>(k, k);
  const arma::mat diag_tt = arma::eye<arma::mat>(tt, tt);
  const arma::mat diag_ktt = arma::eye<arma::mat>(k * tt, k * tt);
  const bool sv = Rcpp::as<bool>(model["sv"]);
  const bool structural = Rcpp::as<bool>(model["structural"]);
  int n_a0 = 0;
  if (structural) {
    n_a0 = k * (k - 1) / 2;
  }
  
  bool covar = false;
  bool varsel = false; // Variable selection
  bool psi_varsel = false;
  bool ssvs = false;
  bool psi_ssvs = false;
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
  
  // Priors & initial values ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  
  Rcpp::List initial = object["initial"];
  Rcpp::List init_coeffs;
  if (use_a) {
    init_coeffs = initial["coefficients"];
  }
  Rcpp::List init_sigma = initial["sigma"];
  
  // Priors - Coefficients
  Rcpp::List priors_coefficients;
  Rcpp::CharacterVector prcoeff_names;
  arma::mat prior_a_mu;
  arma::mat prior_a_Vi;
  Rcpp::List a_prior_varsel;
  arma::vec a_prior_incl, a_tau0, a_tau1, a_tau0sq, a_tau1sq;
  arma::vec a_bvs_lprior_0, a_bvs_lprior_1;
  int a_varsel_n, a_varsel_pos;
  double a_bayes, a_bayes_rand, a_l0, a_l1, a_lambda_draw;
  arma::vec a_post_incl, a_post_mu, a_randu, a_theta0_res, a_theta1_res, a_u0, a_u1, a_varsel_include, a_varsel_include_draw;
  arma::mat a, a_AG, a_lambda, a_post_v, a_theta0, a_theta1, z_bvs;
  if (n_tot > 0) {
    // Priors - Coefficients
    priors_coefficients = priors["coefficients"];
    prcoeff_names = priors_coefficients.names();
    prior_a_mu = Rcpp::as<arma::mat>(priors_coefficients["mu"]);
    prior_a_Vi = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
    
    // Priors - Coefficients - Variable Selection
    if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "ssvs") != prcoeff_names.end()) {
      ssvs = true;
      a_prior_varsel = priors_coefficients["ssvs"];
      a_prior_incl = Rcpp::as<arma::mat>(a_prior_varsel["inprior"]);
      a_tau0 = Rcpp::as<arma::vec>(a_prior_varsel["tau0"]);
      a_tau1 = Rcpp::as<arma::vec>(a_prior_varsel["tau1"]);
      a_tau0sq = arma::square(a_tau0);
      a_tau1sq = arma::square(a_tau1);
      
      if (sv && ssvs) {
        Rcpp::stop("Not allowed to use SSVS with stochastic volatility.");
      }
    }  
    
    if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "bvs") != prcoeff_names.end()) {
      bvs = true;
      a_prior_varsel = priors_coefficients["bvs"];
      a_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
      a_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(a_prior_varsel["inprior"]));
    }
    
    varsel = ssvs || bvs;  
    
    a = Rcpp::as<arma::mat>(init_coeffs["draw"]);
    if (varsel) {
      a_varsel_include = Rcpp::as<arma::vec>(a_prior_varsel["include"]) - 1;
      a_varsel_n = size(a_varsel_include)(0);
      if (ssvs) {
        a_lambda = arma::ones<arma::mat>(n_tot, 1);
      }
      if (bvs) {
        a_lambda = arma::eye<arma::mat>(n_tot, n_tot);
        a_l0 = 0;
        a_l1 = 0;
        a_bayes = 0;
        a_bayes_rand = 0;
        z_bvs = z;
      }
    } 
  }
  
  // Priors - Covar coefficients
  Rcpp::List psi_priors, psi_prior_varsel;
  Rcpp::CharacterVector psi_priors_names;
  arma::vec psi_prior_incl, psi_tau0, psi_tau1, psi_tau0sq, psi_tau1sq, psi_bvs_lprior_0, psi_bvs_lprior_1;
  arma::mat psi_prior_mu, psi_prior_vi;
  int psi_varsel_n, psi_varsel_pos;
  double psi_bayes, psi_bayes_rand, psi_l0, psi_l1, psi_lambda_draw;
  arma::vec psi_post_incl, psi_post_mu, psi_randu, psi_theta0_res, psi_theta1_res, psi_varsel_include, psi_varsel_include_draw, psi_u0, psi_u1, psi_y;
  arma::mat diag_omega_i, diag_covar_omega_i, diag_Psi, psi, Psi, psi_AG, psi_lambda, psi_post_v, psi_theta0, psi_theta1, psi_z, psi_z_bvs;
  if (std::find(priors_names.begin(), priors_names.end(), "psi") != priors_names.end()) {
    covar = true;
    psi_priors = priors["psi"];
    psi_priors_names = psi_priors.names();
    psi_prior_mu = Rcpp::as<arma::mat>(psi_priors["mu"]);
    psi_prior_vi = Rcpp::as<arma::mat>(psi_priors["v_i"]);
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "ssvs") != psi_priors_names.end()) {
      psi_ssvs = true;
      psi_prior_varsel = psi_priors["ssvs"];
      psi_prior_incl = Rcpp::as<arma::mat>(psi_prior_varsel["inprior"]);
      psi_tau0 = Rcpp::as<arma::vec>(psi_prior_varsel["tau0"]);
      psi_tau1 = Rcpp::as<arma::vec>(psi_prior_varsel["tau1"]);
      psi_tau0sq = arma::square(psi_tau0);
      psi_tau1sq = arma::square(psi_tau1);
    }
    if (sv && psi_ssvs) {
      Rcpp::stop("Not allowed to use SSVS with stochastic volatility.");
    }
    
    if (std::find(psi_priors_names.begin(), psi_priors_names.end(), "bvs") != psi_priors_names.end()) {
      psi_bvs = true;
      psi_prior_varsel = psi_priors["bvs"];
      psi_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
      psi_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(psi_prior_varsel["inprior"]));
    }
    
    psi_varsel = psi_ssvs || psi_bvs;
    
    n_psi = k * (k - 1) / 2;
    Psi = arma::eye<arma::mat>(k, k);
    psi_z = arma::zeros<arma::mat>((k - 1) * tt, n_psi);
    if (psi_varsel) {
      psi_varsel_include = Rcpp::as<arma::vec>(psi_prior_varsel["include"]) - 1;
      psi_varsel_n = size(psi_varsel_include)(0);
      if (psi_ssvs) {
        psi_lambda = arma::ones<arma::mat>(n_psi, 1);
      }
      if (psi_bvs) {
        psi_lambda = arma::eye<arma::mat>(n_psi, n_psi);
        psi_l0 = 0;
        psi_l1 = 0;
        psi_bayes = 0;
        psi_bayes_rand = 0;
      }
    }
    diag_covar_omega_i = arma::zeros<arma::mat>(tt * (k - 1), tt * (k - 1));
  }
  
  
  // Priors - Errors
  Rcpp::List sigma_pr = priors["sigma"];
  Rcpp::CharacterVector sigma_names = sigma_pr.names();
  double sigma_post_df;
  arma::vec sigma_post_shape, sigma_prior_rate, sigma_prior_mu;
  arma::mat sigma_prior_scale, sigma_prior_vi;
  bool use_gamma = false;
  arma::vec h_init, sigma_h, u_vec, sigma_post_scale;
  arma::mat h_init_post_v, sigma_h_i, diag_sigma_i_temp;
  arma::vec h_init_post_mu;
  arma::mat sigma_i, h, h_lag, sse, omega_i;
  arma::mat u = y * 0;
  arma::mat diag_sigma_i = arma::zeros(k * tt, k * tt);
  if (sv) {
    sigma_prior_mu = Rcpp::as<arma::vec>(sigma_pr["mu"]);
    sigma_prior_vi = Rcpp::as<arma::mat>(sigma_pr["v_i"]);
    sigma_post_shape = Rcpp::as<arma::vec>(sigma_pr["shape"]) + 0.5 * tt;
    sigma_prior_rate = Rcpp::as<arma::vec>(sigma_pr["rate"]);
    
    h = Rcpp::as<arma::mat>(init_sigma["h"]);
    h_lag = h * 0;
    sigma_h = Rcpp::as<arma::vec>(init_sigma["sigma_h"]);
    h_init = arma::vectorise(h.row(0));
    sigma_i = arma::diagmat(1 / exp(h_init));
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
    
    omega_i = Rcpp::as<arma::mat>(init_sigma["sigma_i"]);
    sigma_i = omega_i;
  }
  diag_sigma_i.diag() = arma::repmat(sigma_i.diag(), tt, 1);
  if (covar | sv) {
    diag_omega_i = diag_sigma_i;
  }

  // Storage objects
  const int iter = Rcpp::as<int>(model["iterations"]);
  const int burnin = Rcpp::as<int>(model["burnin"]);
  const int draws = iter + burnin;
  int pos_draw;
  const int a_pos_start = 0;
  const int a_pos_end = n_a - 1;
  const int b_pos_start = n_a;
  const int b_pos_end = n_a + n_b - 1;
  const int c_pos_start = n_a + n_b;
  const int c_pos_end = n_a + n_b + n_c - 1;
  const int a0_pos_start = n_a + n_b + n_c;
  const int a0_pos_end = n_a + n_b + n_c + n_a0 - 1;

  arma::mat draws_a0 = arma::zeros<arma::mat>(n_a0, iter);
  arma::mat draws_a = arma::zeros<arma::mat>(n_a, iter);
  arma::mat draws_b = arma::zeros<arma::mat>(n_b, iter);
  arma::mat draws_c = arma::zeros<arma::mat>(n_c, iter);
  arma::mat draws_sigma, draws_sigma_sigma;
  if (sv) {
    draws_sigma = arma::zeros<arma::mat>(k * k * tt, iter);
    draws_sigma_sigma = arma::zeros<arma::mat>(k * k, iter);
  } else {
    draws_sigma = arma::zeros<arma::mat>(k * k, iter);
  }
  
  arma::vec lambda_vec, psi_lambda_vec;
  arma::mat draws_lambda_a0, draws_lambda_a, draws_lambda_b, draws_lambda_c;
  if (varsel) {
    if (structural) {
      draws_lambda_a0 = arma::zeros<arma::mat>(n_a0, iter); 
    }
    draws_lambda_a = arma::zeros<arma::mat>(n_a, iter);
    draws_lambda_b = arma::zeros<arma::mat>(n_b, iter);
    draws_lambda_c = arma::zeros<arma::mat>(n_c, iter);
  }
  if (covar && psi_varsel) {
    draws_lambda_a0 = arma::zeros<arma::mat>(n_psi, iter);  
  }

  // Start Gibbs sampler
  for (int draw = 0; draw < draws; draw++) {
    
    if (draw % 20 == 0) { // Check for user interruption every now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Draw coefficients ----
    if (use_a) {
      
      if (bvs) {
        z = z_bvs * a_lambda;
      }
      a_post_v = prior_a_Vi + arma::trans(z) * diag_sigma_i * z;
      a_post_mu = arma::solve(a_post_v, prior_a_Vi * prior_a_mu + arma::trans(z) * diag_sigma_i * yvec);
      a = a_post_mu + arma::solve(arma::chol(a_post_v), arma::randn(n_tot));
      
      if (varsel) {
        
        // Reorder positions of variable selection
        a_varsel_include_draw = shuffle(a_varsel_include);
        
        if (ssvs) {
          // Obtain inclusion posterior
          a_u0 = 1 / a_tau0 % arma::exp(-(arma::square(a) / (2 * a_tau0sq))) % (1 - a_prior_incl);
          a_u1 = 1 / a_tau1 % arma::exp(-(arma::square(a) / (2 * a_tau1sq))) % a_prior_incl;
          a_post_incl = a_u1 / (a_u0 + a_u1);
          
          // Draw inclusion parameters in random order
          for (int i = 0; i < a_varsel_n; i++){
            a_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, a_post_incl(a_varsel_include_draw(i))));
            a_lambda(a_varsel_include_draw(i), 0) = a_lambda_draw;
            if (a_lambda_draw == 0) {
              prior_a_Vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau0sq(a_varsel_include_draw(i));
            } else {
              prior_a_Vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau1sq(a_varsel_include_draw(i));
            }
          }
          lambda_vec = arma::vectorise(a_lambda);
        }
        
        if (bvs) {
          z = z_bvs;
          a_AG = a_lambda * a;
          for (int j = 0; j < a_varsel_n; j++){
            a_varsel_pos = a_varsel_include_draw(j);
            a_randu = arma::log(arma::randu<arma::vec>(1));
            if (a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) >= a_bvs_lprior_1(a_varsel_pos)){continue;}
            if (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) >= a_bvs_lprior_0(a_varsel_pos)){continue;}
            if ((a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) < a_bvs_lprior_1(a_varsel_pos)) || (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) < a_bvs_lprior_0(a_varsel_pos))){
              a_theta0 = a_AG;
              a_theta1 = a_AG;
              a_theta0.row(a_varsel_pos) = 0;
              a_theta1.row(a_varsel_pos) = a.row(a_varsel_pos);
              a_theta0_res = yvec - z * a_theta0;
              a_theta1_res = yvec - z * a_theta1;
              a_l0 = -arma::as_scalar(trans(a_theta0_res) * diag_sigma_i * a_theta0_res) / 2 + arma::as_scalar(a_bvs_lprior_0(a_varsel_pos));
              a_l1 = -arma::as_scalar(trans(a_theta1_res) * diag_sigma_i * a_theta1_res) / 2 + arma::as_scalar(a_bvs_lprior_1(a_varsel_pos));
              a_bayes = a_l1 - a_l0;
              a_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
              if (a_bayes >= a_bayes_rand){
                a_lambda(a_varsel_pos, a_varsel_pos) = 1;
              } else {
                a_lambda(a_varsel_pos, a_varsel_pos) = 0;
              }
            }
          }
          a = a_lambda * a;
          lambda_vec = a_lambda.diag();
        }
      }
      
      u_vec = yvec - z * a;
      u = arma::reshape(u_vec, k, tt);  
    } else {
      u = y;
    }
    
    // Covariances
    if (covar) {
      
      // Prepare data
      psi_y = arma::vectorise(u.rows(1, k - 1));
      for (int i = 1; i < k; i++) {
        for (int j = 0; j < tt; j++) {
          psi_z.submat(j * (k - 1) + i - 1,
                       i * (i - 1) / 2,
                       j * (k - 1) + i - 1,
                       (i + 1) * i / 2 - 1) = -arma::trans(u.submat(0, j, i - 1, j));
          
          diag_covar_omega_i(j * (k - 1) + i - 1, j * (k - 1) + i - 1) = diag_omega_i(j * k + i, j * k + i);
        }
      }
      
      if (psi_bvs) {
        psi_z_bvs = psi_z;
        psi_z = psi_z_bvs * psi_lambda;
      }
      psi_post_v = psi_prior_vi + arma::trans(psi_z) * diag_covar_omega_i * psi_z;
      psi_post_mu = arma::solve(psi_post_v, psi_prior_vi * psi_prior_mu + arma::trans(psi_z) * diag_covar_omega_i * psi_y);
      psi = psi_post_mu + arma::solve(arma::chol(psi_post_v), arma::randn(n_psi));
      
      if (psi_varsel) {
        
        // Reorder positions of variable selection
        psi_varsel_include_draw = shuffle(psi_varsel_include);
        
        if (psi_ssvs) {
          // Obtain inclusion posterior
          psi_u0 = 1 / psi_tau0 % arma::exp(-(arma::square(psi) / (2 * psi_tau0sq))) % (1 - psi_prior_incl);
          psi_u1 = 1 / psi_tau1 % arma::exp(-(arma::square(psi) / (2 * psi_tau1sq))) % psi_prior_incl;
          psi_post_incl = psi_u1 / (psi_u0 + psi_u1);
          
          // Draw inclusion parameters in random order
          for (int i = 0; i < psi_varsel_n; i++){
            psi_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, psi_post_incl(psi_varsel_include_draw(i))));
            psi_lambda(psi_varsel_include_draw(i), 0) = psi_lambda_draw;
            if (psi_lambda_draw == 0) {
              psi_prior_vi(psi_varsel_include_draw(i), psi_varsel_include_draw(i)) = 1 / psi_tau0sq(psi_varsel_include_draw(i));
            } else {
              psi_prior_vi(psi_varsel_include_draw(i), psi_varsel_include_draw(i)) = 1 / psi_tau1sq(psi_varsel_include_draw(i));
            }
          }
          psi_lambda_vec = arma::vectorise(psi_lambda);
        }
        
        if (psi_bvs) {
          psi_z = psi_z_bvs;
          psi_AG = psi_lambda * psi;
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
          psi = psi_lambda * psi;
          psi_lambda_vec = psi_lambda.diag();
        }
      }
      
      for (int i = 1; i < k; i++) {
        Psi.submat(i, 0, i, i - 1) = arma::trans(psi.submat(i * (i - 1) / 2, 0, (i + 1) * i / 2 - 1, 0));
      }
      u = Psi * u;
    }
    
    if (sv) {
      
      // Draw variances  
      for (int i = 0; i < k; i++) {
        h.col(i) = bvartools::stoch_vol(u.row(i).t(), h.col(i), sigma_h(i), h_init(i));
      }
      diag_omega_i.diag() = 1 / exp(arma::vectorise(h.t()));
      if (covar) {
        diag_Psi = arma::kron(diag_tt, Psi);
        diag_sigma_i = arma::trans(diag_Psi) * diag_omega_i * diag_Psi;
      } else {
        diag_sigma_i = diag_omega_i;
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
        if (covar) {
          diag_omega_i = arma::kron(diag_tt, omega_i);
          sigma_i = arma::trans(Psi) * omega_i * Psi;
        } else {
          sigma_i = omega_i;
        }
        
      } else {
        sigma_i = arma::wishrnd(arma::solve(sigma_prior_scale + u * u.t(), diag_k), sigma_post_df);
      }
      
      diag_sigma_i = arma::kron(diag_tt, sigma_i);
      
    }
    
    // Store draws
    if (draw >= burnin) {
      
      pos_draw = draw - burnin;
      
      if (sv) {
        for (int i = 0; i < tt; i ++) {
          draws_sigma.submat(i * n_sigma, pos_draw, (i + 1) * n_sigma - 1, pos_draw) = arma::vectorise(arma::solve(diag_sigma_i.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1), diag_k));
        }
        draws_sigma_sigma.col(pos_draw) = arma::vectorise(arma::diagmat(sigma_h));
      } else {
        draws_sigma.col(pos_draw) = arma::vectorise(arma::solve(sigma_i, diag_k));
      }
      
      if (psi_varsel) {
        draws_lambda_a0.col(pos_draw) = psi_lambda_vec;
      }
      
      if (n_a > 0) {
        draws_a.col(pos_draw) = arma::vectorise(a.rows(a_pos_start, a_pos_end));
        if (varsel) {
          draws_lambda_a.col(pos_draw) = lambda_vec.subvec(a_pos_start, a_pos_end);
        }
      }
      if (n_b > 0) {
        draws_b.col(pos_draw) = arma::vectorise(a.rows(b_pos_start, b_pos_end));
        if (varsel) {
          draws_lambda_b.col(pos_draw) = lambda_vec.subvec(b_pos_start, b_pos_end);
        }
      }
      if (n_c > 0) {
        draws_c.col(pos_draw) = arma::vectorise(a.rows(c_pos_start, c_pos_end));
        if (varsel) {
          draws_lambda_c.col(pos_draw) = lambda_vec.subvec(c_pos_start, c_pos_end);
        }
      }
      if (structural) {
        draws_a0.col(pos_draw) = arma::vectorise(a.rows(a0_pos_start, a0_pos_end));
        if (varsel) {
          draws_lambda_a0.col(pos_draw) = lambda_vec.subvec(a0_pos_start, a0_pos_end);
        }
      }
    }
  } // End loop

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("b") = R_NilValue,
                                             Rcpp::Named("c") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);

  if (n_a > 0) {
    if (varsel) {
      posteriors["a"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a,
                                                     Rcpp::Named("lambda") = draws_lambda_a));
    } else {
      posteriors["a"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a));
    }
  }

  if (n_b > 0) {
    if (varsel) {
      posteriors["b"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_b,
                                                     Rcpp::Named("lambda") = draws_lambda_b));
    } else {
      posteriors["b"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_b));
    }
  }

  if (n_c > 0) {
    if (varsel) {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c,
                                                      Rcpp::Named("lambda") = draws_lambda_c));
    } else {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_c));
    }
  }

  if (structural) {
    if (varsel) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0,
                                                       Rcpp::Named("lambda") = draws_lambda_a0));
    } else {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a0));
    }
  }

  if (psi_varsel) {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("lambda") = draws_lambda_a0));
    }
  } else {
    if (sv) {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                     Rcpp::Named("sigma") = draws_sigma_sigma)); 
    } else {
      posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma));
    }
  }

  Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = object["data"],
                                         Rcpp::Named("model") = object["model"],
                                         Rcpp::Named("priors") = object["priors"],
                                         Rcpp::Named("posteriors") = posteriors);

  return result;
  // return Rcpp::List::create(Rcpp::Named("test") = sigma_i);
}