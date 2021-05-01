#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(.bvecalg)]]
Rcpp::List bvecalg(Rcpp::List object) {
  
  // Initialise variables
  Rcpp::List data = object["data"];
  arma::mat y = arma::trans(Rcpp::as<arma::mat>(data["Y"]));
  arma::mat w = arma::trans(Rcpp::as<arma::mat>(data["W"]));
  Rcpp::Nullable<Rcpp::List> x_test = data["X"];
  arma::mat x;
  if (x_test.isNotNull()) {
    x = arma::trans(Rcpp::as<arma::mat>(data["X"])); 
  }
  arma::mat yvec = arma::vectorise(y);
  Rcpp::List model = object["model"];
  Rcpp::CharacterVector model_names = model.names();
  Rcpp::List endogen = model["endogen"];
  Rcpp::CharacterVector endo_names = Rcpp::as<Rcpp::CharacterVector>(endogen["variables"]);
  
  // Define useful variables
  int tt = y.n_cols;
  int k = y.n_rows;
  int p = Rcpp::as<int>(endogen["lags"]);
  int n_gamma = 0;
  if (p > 1) {
    n_gamma = k * k * (p - 1); 
  }
  int m = 0;
  int s = 0;
  int n_upsilon = 0;
  int n_x = 0;
  if (x_test.isNotNull()) {
    n_x = x.n_rows; 
  }
  
  bool use_rr = false;
  arma::mat alpha, Alpha, beta, Beta, pi, y_beta;
  int r = Rcpp::as<int>(model["rank"]);
  int n_w = 0;
  if (r > 0) {
    use_rr = true;
    n_w = w.n_rows;
    alpha = arma::zeros<arma::mat>(k, r);
    beta = arma::zeros<arma::mat>(n_w, r);
    beta.submat(0, 0, r - 1, r - 1) = arma::eye<arma::mat>(r, r);
    pi = arma::zeros<arma::mat>(k, n_w);
  }
  int n_alpha = k * r;
  int n_beta = n_w * r;
  
  int n_a0 = 0;
  bool structural = Rcpp::as<bool>(model["structural"]);
  if (structural && k > 1) {
    n_a0 = k * (k - 1) / 2;
  }
  
  int n_tot = n_alpha + k * n_x + n_a0;
  arma::mat diag_k = arma::eye<arma::mat>(k, k);
  arma::mat diag_beta;
  if (use_rr) {
    diag_beta = arma::eye<arma::mat>(n_beta, n_beta);
  }
  arma::mat diag_tot = arma::eye<arma::mat>(n_tot, n_tot);
  arma::mat diag_tt = arma::eye<arma::mat>(tt, tt);
  
  arma::mat z = arma::zeros<arma::mat>(tt * k, n_alpha + k * n_x + n_a0);
  arma::mat BB_sqrt, z_beta, svd_U, svd_V;
  arma::vec svd_s;
  if (use_rr) {
    z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
    if (n_x > 0) {
      z.cols(n_alpha, n_alpha + n_x * k + n_a0 - 1) = Rcpp::as<arma::mat>(data["SUR"]).cols(k * n_w, k * (n_w + n_x) + n_a0 - 1);
    }
    z_beta = arma::zeros<arma::mat>(k * tt, n_beta);
  } else {
    z = Rcpp::as<arma::mat>(data["SUR"]).cols(k * n_w, k * (n_w + n_x) + n_a0 - 1);
  }
  
  bool varsel = false;
  bool ssvs = false;
  bool bvs = false;
  bool vs_covar = false;
  
  // OLS estimates for initial values
  arma::mat gamma = arma::zeros<arma::mat>(n_tot, 1);
  arma::mat u = arma::reshape(yvec - z * gamma, k, tt);
  arma::mat sigma = u * u.t() / tt;
  arma::mat sigma_i = arma::inv(sigma);
  arma::mat diag_sigma_i;
  
  // Exogenous variables
  Rcpp::CharacterVector exogen_names;
  Rcpp::List exogen;
  if (std::find(model_names.begin(), model_names.end(), "exogen") != model_names.end()) {
    exogen = model["exogen"];
    exogen_names = Rcpp::as<Rcpp::CharacterVector>(exogen["variables"]);
    m = exogen_names.length();
    s = Rcpp::as<int>(exogen["lags"]);
    n_upsilon = k * m * s;
  }
  
  // Deterministic terms
  int n_r = 0;
  int n_ur = 0;
  int n_c_ur = 0;
  Rcpp::CharacterVector determ_names, det_names_r, det_names_ur;
  Rcpp::List determ;
  if (std::find(model_names.begin(), model_names.end(), "deterministic") != model_names.end()) {
    determ = model["deterministic"];
    Rcpp::CharacterVector determ_names = determ.names();
    if (std::find(determ_names.begin(), determ_names.end(), "restricted") != determ_names.end()) {
      det_names_r = Rcpp::as<Rcpp::CharacterVector>(determ["restricted"]);
      n_r = det_names_r.length();
    }
    if (std::find(determ_names.begin(), determ_names.end(), "unrestricted") != determ_names.end()) {
      det_names_ur = Rcpp::as<Rcpp::CharacterVector>(determ["unrestricted"]);
      n_ur = det_names_ur.length();
      n_c_ur = k * n_ur;
    }
  }
  
  // Priors ----
  Rcpp::List priors = object["priors"];
  Rcpp::CharacterVector priors_names = priors.names();
  // Cointegration
  Rcpp::List priors_cointegration;
  Rcpp::CharacterVector prcoint_names;
  if (use_rr) {
    priors_cointegration = priors["cointegration"];
    prcoint_names = priors_cointegration.names();
  }
  // Non-cointegration
  Rcpp::List priors_coefficients;
  Rcpp::CharacterVector prcoeff_names;
  if (std::find(priors_names.begin(), priors_names.end(), "noncointegration") != priors_names.end()) { 
    priors_coefficients = priors["noncointegration"];
    prcoeff_names = priors_coefficients.names(); 
  } else {
    if (!use_rr) {
      Rcpp::stop("Cointegration rank is zero and no non-cointegration priors provided.");
    }
  }
  // Errors
  Rcpp::List sigma_pr = priors["sigma"];
  Rcpp::CharacterVector sigma_names = sigma_pr.names();
  
  // Priors - Coefficients
  double coint_v_i;
  arma::mat prior_a_mu, prior_a_Vi, post_a_mu, post_a_Vi, post_beta_mu, post_beta_Vi, p_tau_i;
  if (use_rr) {
    prior_a_mu = arma::zeros<arma::mat>(n_tot, 1);
    prior_a_Vi = arma::zeros<arma::mat>(n_tot, n_tot);
    coint_v_i = Rcpp::as<double>(priors_cointegration["v_i"]);
    p_tau_i = Rcpp::as<arma::mat>(priors_cointegration["p_tau_i"]);
    if (n_x > 0) {
      prior_a_mu.submat(n_alpha, 0, n_tot - 1, 0) = Rcpp::as<arma::mat>(priors_coefficients["mu"]); 
      prior_a_Vi.submat(n_alpha, n_alpha, n_tot - 1, n_tot - 1) = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
    }
  } else {
    prior_a_mu = Rcpp::as<arma::mat>(priors_coefficients["mu"]);
    prior_a_Vi = Rcpp::as<arma::mat>(priors_coefficients["v_i"]);
  }
  post_a_mu = prior_a_mu * 0;
  post_a_Vi = prior_a_Vi * 0;
  
  // Priors - Coefficient - Variables selection
  arma::vec a_prior_incl, a_varsel_include, a_varsel_include_draw;
  arma::mat a_lambda;
  int a_varsel_n, a_varsel_pos;
  
  // Priors - Coefficient - SSVS
  Rcpp::List a_pr_ssvs;
  arma::vec post_a_incl, a_tau0, a_tau1, a_tau0sq, a_tau1sq, a_u0, a_u1;
  double a_lambda_draw;
  if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "ssvs") != prcoeff_names.end()) {
    if (n_tot - n_alpha > 0) {
      ssvs = true;
      a_pr_ssvs = priors_coefficients["ssvs"];
      a_prior_incl = arma::zeros<arma::mat>(n_tot, 1);
      a_prior_incl.submat(n_alpha, 0, n_tot - 1, 0) = Rcpp::as<arma::mat>(a_pr_ssvs["inprior"]);
      a_tau0 = arma::zeros<arma::vec>(n_tot);
      a_tau1 = arma::zeros<arma::vec>(n_tot);
      a_tau0.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(a_pr_ssvs["tau0"]);
      a_tau1.subvec(n_alpha, n_tot - 1) = Rcpp::as<arma::vec>(a_pr_ssvs["tau1"]);
      a_tau0sq = arma::square(a_tau0);
      a_tau1sq = arma::square(a_tau1);
      a_varsel_include = Rcpp::as<arma::vec>(a_pr_ssvs["include"]) - 1 + n_alpha;
      a_varsel_n = size(a_varsel_include)(0);
      a_lambda = arma::ones<arma::mat>(n_tot, 1);
    } else {
      Rcpp::stop("Model does not contain non-cointegration variables for BVS.");
    }
  }
  
  
  // Priors - Coefficient - BVS
  Rcpp::List a_pr_bvs;
  arma::mat a_AG, a_theta0, a_theta1, z_bvs;
  arma::vec a_bvs_lprior_0, a_bvs_lprior_1, a_bvs_l0_res, a_bvs_l1_res, a_randu;
  double a_l0, a_l1, a_bayes, a_bayes_rand;
  
  if (std::find(prcoeff_names.begin(), prcoeff_names.end(), "bvs") != prcoeff_names.end()) {
    if (n_tot - n_alpha > 0) {
      bvs = true;
      a_pr_bvs = priors_coefficients["bvs"];
      a_bvs_lprior_0 = arma::zeros<arma::vec>(n_tot);
      a_bvs_lprior_1 = arma::zeros<arma::vec>(n_tot);
      a_bvs_lprior_0.submat(n_alpha, 0, n_tot - 1, 0) = arma::log(1 - Rcpp::as<arma::vec>(a_pr_bvs["inprior"]));
      a_bvs_lprior_1.submat(n_alpha, 0, n_tot - 1, 0) = arma::log(Rcpp::as<arma::vec>(a_pr_bvs["inprior"]));
      a_varsel_include = Rcpp::as<arma::vec>(a_pr_bvs["include"]) - 1 + n_alpha;
      a_varsel_n = size(a_varsel_include)(0);
      a_lambda = arma::eye<arma::mat>(n_tot, n_tot);
      a_l0 = 0;
      a_l1 = 0;
      a_bayes = 0;
      a_bayes_rand = 0;
      z_bvs = z;
    } else {
      Rcpp::stop("Model does not contain non-cointegration variables for BVS.");
    }
  }
  
  varsel = ssvs || bvs;
  
  // Priors - Errors
  double prior_sigma_df, post_sigma_df, sigma_prior_shape, sigma_post_shape, sigma_post_scale_j;
  arma::mat sigma_prior_rate, sigma_prior_scale, psi_prior_Vi, psi_prior_Vi_j, sse, u_vec;
  bool use_gamma = false;
  if (std::find(sigma_names.begin(), sigma_names.end(), "df") != sigma_names.end()) {
    prior_sigma_df = Rcpp::as<double>(sigma_pr["df"]);
    post_sigma_df = prior_sigma_df + tt;
    sigma_prior_scale = Rcpp::as<arma::mat>(sigma_pr["scale"]);
  }
  if (std::find(sigma_names.begin(), sigma_names.end(), "shape") != sigma_names.end()) {
    use_gamma = true;
    sigma_prior_shape = Rcpp::as<double>(sigma_pr["shape"]);
    sigma_post_shape = sigma_prior_shape + tt / 2;
    sigma_prior_rate = Rcpp::as<arma::mat>(sigma_pr["rate"]);
    sigma_i = arma::diagmat(sigma_i.diag());
  }
  
  // Priors - Errors - Covar
  Rcpp::List psi_pr;
  arma::mat Psi;
  arma::vec psi_u0, psi_u1, psi_prior_incl, psi_post_incl;
  
  // Priors - Errors - SSVS
  int psi_pos1, psi_pos2;
  double psi_lambda_draw;
  arma::mat psi_j, s_j, sse_j, psi_post_V_j, psi_post_mu_j;
  arma::vec psi_lambda, psi_tau0sq, psi_tau0, psi_tau1, psi_tau1sq;
  if (std::find(sigma_names.begin(), sigma_names.end(), "ssvs") != sigma_names.end()) {
    vs_covar = true;
    psi_pr = sigma_pr["ssvs"];
    psi_prior_Vi = Rcpp::as<arma::mat>(sigma_pr["v_i"]);
    psi_prior_incl = Rcpp::as<arma::vec>(psi_pr["inprior"]);
    psi_tau0 = Rcpp::as<arma::vec>(psi_pr["tau0"]);
    psi_tau1 = Rcpp::as<arma::vec>(psi_pr["tau1"]);
    psi_tau0sq = arma::square(psi_tau0);
    psi_tau1sq = arma::square(psi_tau1);
    psi_lambda = arma::ones<arma::vec>(k * (k - 1) / 2);
    psi_u0 = psi_tau0 * 0;
    psi_u1 = psi_tau1 * 0;
    
    Psi = arma::eye<arma::mat>(k, k);
    Psi.diag() = sqrt(sigma_i.diag());
  }
  
  // Priors - Errors - BVS
  arma::mat cov_lambda, diag_omega_i, omega_i, psi, psi_post_mu, psi_post_V, psi_prior_mu, z_covar, z_covar_bvs;
  arma::mat cov_AG, cov_bvs_l0_res, cov_bvs_l1_res, cov_bvs_lprior_0, cov_bvs_lprior_1, cov_theta0, cov_theta1;
  arma::vec cov_randu;
  double cov_l0, cov_l1, cov_bayes, cov_bayes_rand;
  if (std::find(sigma_names.begin(), sigma_names.end(), "bvs") != sigma_names.end()) {
    vs_covar = true;
    n_a0 = k * (k - 1) / 2;
    psi_pr = sigma_pr["bvs"];
    cov_bvs_lprior_0 = arma::log(1 - Rcpp::as<arma::vec>(psi_pr["inprior"]));
    cov_bvs_lprior_1 = arma::log(Rcpp::as<arma::vec>(psi_pr["inprior"]));
    psi_prior_mu = Rcpp::as<arma::mat>(sigma_pr["mu"]);
    psi_prior_Vi = Rcpp::as<arma::mat>(sigma_pr["v_i"]);
    z_covar = arma::zeros<arma::mat>(tt * k, n_a0);
    z_covar_bvs = z_covar;
    cov_lambda = arma::eye<arma::mat>(n_a0, n_a0);
    omega_i = arma::diagmat(sigma_i.diag());
    
    Psi = arma::eye<arma::mat>(k, k);
  }
  
  if (vs_covar && k == 1) {
    Rcpp::stop("Not allowed to apply BVS or SSVS to covariances when there is only one endogenous variable.");
  }
  
  if (vs_covar && structural) {
    Rcpp::stop("Not allowed to apply BVS or SSVS to covariances when a structural model is estimated.");
  }
  
  // Storage objects
  int iter = Rcpp::as<int>(model["iterations"]);
  int burnin = Rcpp::as<int>(model["burnin"]);
  int draws = iter + burnin;
  int pos_draw = 0;
  
  arma::mat draws_a = arma::zeros<arma::mat>(n_tot, iter);
  
  arma::mat draws_alpha;
  arma::mat draws_beta;
  arma::mat draws_pi;
  if (use_rr) {
    draws_alpha = arma::zeros<arma::mat>(n_alpha, iter);
    draws_beta = arma::zeros<arma::mat>(n_beta, iter);
    draws_pi = arma::zeros<arma::mat>(k * n_w, iter);
  }
  arma::mat draws_sigma = arma::zeros<arma::mat>(k * k, iter);
  
  arma::mat draws_a_lambda, draws_cov_lambda;
  if (varsel) {
    draws_a_lambda = arma::zeros<arma::mat>(n_tot, iter);
  }
  if (vs_covar) {
    draws_cov_lambda = arma::zeros<arma::mat>(k * (k - 1) / 2, iter);
  }
  
  // Rcpp::List bla = Rcpp::List::create(Rcpp::Named("test") = pos_draw);
  // return bla;
  
  for (int draw = 0; draw < draws; draw++) {
    
    if (draw % 1000 == 0) { // Check for user interuption ever now and then
      Rcpp::checkUserInterrupt();
    }
    
    // Draw error
    if (use_rr) {
      u_vec = yvec - arma::vectorise(pi * w) - z.cols(n_alpha, n_tot - 1) * gamma.rows(n_alpha, n_tot - 1); 
    } else {
      u_vec = yvec - z * gamma;
    }
    u = arma::reshape(u_vec, k, tt);
    if (use_gamma) {
      
      sse = u * u.t();
      
      if (vs_covar) {
        
        if (ssvs) {
          for (int i = 0; i < k; i++) {
            
            psi_pos1 = i * (i - 1) / 2;
            psi_pos2 = (i + 1) * i / 2 - 1;
            
            if (i > 0) {
              psi_prior_Vi_j = psi_prior_Vi.submat(psi_pos1, psi_pos1, psi_pos2, psi_pos2);
              sse_j = sse.submat(0, 0, i - 1, i - 1);
              s_j = sse.submat(0, i, i - 1, i);
              psi_post_V_j = arma::inv(psi_prior_Vi_j + sse_j);
              psi_post_mu_j = -arma::as_scalar(Psi(i, i)) * psi_post_V_j * s_j;
              psi_j = arma::mvnrnd(psi_post_mu_j, psi_post_V_j);
              Psi.submat(0, i, i - 1, i) = psi_j;
              
              // Obtain inclusion posterior
              psi_u0.subvec(psi_pos1, psi_pos2) = 1 / psi_tau0.subvec(psi_pos1, psi_pos2) % arma::exp(-(arma::square(psi_j) / (2 * psi_tau0sq.subvec(psi_pos1, psi_pos2)))) % (1 - psi_prior_incl.subvec(psi_pos1, psi_pos2));
              psi_u1.subvec(psi_pos1, psi_pos2) = 1 / psi_tau1.subvec(psi_pos1, psi_pos2) % arma::exp(-(arma::square(psi_j) / (2 * psi_tau1sq.subvec(psi_pos1, psi_pos2)))) % psi_prior_incl.subvec(psi_pos1, psi_pos2);
              psi_post_incl = psi_u1.subvec(psi_pos1, psi_pos2) / (psi_u0.subvec(psi_pos1, psi_pos2) + psi_u1.subvec(psi_pos1, psi_pos2));
              
              // Draw inclusion parameters in random order
              for (int j = 0; j < i; j++){
                psi_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, psi_post_incl(j)));
                psi_lambda(psi_pos1 + j) = psi_lambda_draw;
                if (psi_lambda_draw == 0) {
                  psi_prior_Vi(psi_pos1 + j, psi_pos1 + j) = 1 / psi_tau0sq(psi_pos1 + j);
                } else {
                  psi_prior_Vi(psi_pos1 + j, psi_pos1 + j) = 1 / psi_tau1sq(psi_pos1 + j);
                }
              }
            }
            // Draw variances
            if (i == 0) {
              sigma_post_scale_j = 1 / arma::as_scalar(sigma_prior_rate(0, 0) + sse(0, 0) / 2);
            } else {
              sigma_post_scale_j = 1 / arma::as_scalar(sigma_prior_rate(i, i) + (sse(i, i) - s_j.t() * psi_post_V_j * s_j) / 2);
            }
            
            Psi(i, i) =  sqrt(arma::randg<double>(arma::distr_param(sigma_post_shape, sigma_post_scale_j)));
          }
          sigma_i = Psi * Psi.t();
        }
        
        if (bvs) {
          
          for (int j = 0; j < tt; j++) {
            for (int l = 1; l < k; l++) {
              z_covar_bvs.submat(j * k + l, l * (l - 1) / 2, j * k + l, (l + 1) * l / 2 - 1) = -arma::trans(u.submat(0, j, l - 1, j));
            }
          }
          z_covar = z_covar_bvs * cov_lambda;
          diag_omega_i = arma::kron(diag_tt, omega_i);
          psi_post_V = arma::inv(psi_prior_Vi + arma::trans(z_covar) * diag_omega_i *  z_covar);
          psi_post_mu = psi_post_V * (psi_prior_Vi * psi_prior_mu + arma::trans(z_covar) * diag_omega_i * u_vec);
          psi = arma::mvnrnd(psi_post_mu, psi_post_V);
          
          // Add BVS here
          z_covar = z_covar_bvs;
          cov_AG = cov_lambda * psi;
          for (int j = 0; j < n_a0; j++){
            cov_randu = arma::log(arma::randu<arma::vec>(1));
            if (cov_lambda(j, j) == 1 && cov_randu(0) >= cov_bvs_lprior_1(j)){continue;}
            if (cov_lambda(j, j) == 0 && cov_randu(0) >= cov_bvs_lprior_0(j)){continue;}
            if ((cov_lambda(j, j) == 1 && cov_randu(0) < cov_bvs_lprior_1(j)) || (cov_lambda(j, j) == 0 && cov_randu(0) < cov_bvs_lprior_0(j))){
              cov_theta0 = cov_AG;
              cov_theta1 = cov_AG;
              cov_theta1.row(j) = psi.row(j);
              cov_theta0.row(j) = 0;
              cov_bvs_l0_res = u_vec - z_covar * cov_theta0;
              cov_bvs_l1_res = u_vec - z_covar * cov_theta1;
              cov_l0 = -arma::as_scalar(trans(cov_bvs_l0_res) * diag_omega_i * cov_bvs_l0_res) / 2 + arma::as_scalar(cov_bvs_lprior_0(j));
              cov_l1 = -arma::as_scalar(trans(cov_bvs_l1_res) * diag_omega_i * cov_bvs_l1_res) / 2 + arma::as_scalar(cov_bvs_lprior_1(j));
              cov_bayes = cov_l1 - cov_l0;
              cov_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
              if (cov_bayes >= cov_bayes_rand){
                cov_lambda(j, j) = 1;
              } else {
                cov_lambda(j, j) = 0;
              }
            }
          }
          psi = cov_lambda * psi;
          
          for (int j = 1; j < k; j++) {
            Psi.submat(j, 0, j, j - 1) = arma::trans(psi.submat(j * (j - 1) / 2, 0, (j + 1) * j / 2 - 1, 0));
          }
          u = Psi * u;
          sse = u * u.t();
          
          // Draw variances
          for (int j = 0; j < k; j++) {
            omega_i(j, j) = arma::randg<double>(arma::distr_param(sigma_post_shape, 1 / arma::as_scalar(sigma_prior_rate(j, j) + sse(j, j) / 2)));
          }
          
          sigma_i = arma::trans(Psi) * omega_i * Psi;
        }
        
      } else {
        for (int i = 0; i < k; i++) {
          sigma_i(i, i) = arma::randg<double>(arma::distr_param(sigma_post_shape, 1 / arma::as_scalar(sigma_prior_rate(i, i) + sse(i, i) / 2)));
        } 
      }
    } else {
      // Wishart
      if (use_rr) {
        sigma_i = arma::wishrnd(arma::solve(coint_v_i * alpha * (beta.t() * p_tau_i * beta) * alpha.t() + u * u.t(), diag_k), post_sigma_df); 
      } else {
        sigma_i = arma::wishrnd(arma::solve(sigma_prior_scale + u * u.t(), diag_k), post_sigma_df);
      }
    }
    
    // Draw coefficients ----
    diag_sigma_i = arma::kron(diag_tt, sigma_i);
    if (use_rr) { // Update priors for alpha
      post_a_Vi.submat(0, 0, n_alpha - 1, n_alpha - 1) = arma::kron(coint_v_i * (arma::trans(beta) * p_tau_i * beta), sigma_i);
      if (bvs) {
        z_bvs.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
      } else {
        z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta) * w), diag_k);
      }
    }
    if (bvs) {
      z = z_bvs * a_lambda;
    }
    post_a_Vi = arma::solve(prior_a_Vi + arma::trans(z) * diag_sigma_i * z, diag_tot);
    post_a_mu = post_a_Vi * (prior_a_Vi * prior_a_mu + arma::trans(z) * diag_sigma_i * yvec);
    gamma = arma::mvnrnd(post_a_mu, post_a_Vi);
    
    if (varsel) {
      a_varsel_include_draw = shuffle(a_varsel_include); // Reorder positions of variable selection
      
      if (ssvs) {
        // Obtain inclusion posterior
        a_u0 = 1 / a_tau0 % arma::exp(-(arma::square(gamma) / (2 * a_tau0sq))) % (1 - a_prior_incl);
        a_u1 = 1 / a_tau1 % arma::exp(-(arma::square(gamma) / (2 * a_tau1sq))) % a_prior_incl;
        post_a_incl = a_u1 / (a_u0 + a_u1);
        
        // Draw inclusion parameters in random order
        for (int i = 0; i < a_varsel_n; i++){
          a_lambda_draw = Rcpp::as<double>(Rcpp::rbinom(1, 1, post_a_incl(a_varsel_include_draw(i))));
          a_lambda(a_varsel_include_draw(i), 0) = a_lambda_draw;
          if (a_lambda_draw == 0) {
            prior_a_Vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau0sq(a_varsel_include_draw(i));
          } else {
            prior_a_Vi(a_varsel_include_draw(i), a_varsel_include_draw(i)) = 1 / a_tau1sq(a_varsel_include_draw(i));
          }
        }
      }
      
      if (bvs) {
        z = z_bvs;
        a_AG = a_lambda * gamma;
        for (int j = 0; j < a_varsel_n; j++){
          a_varsel_pos = a_varsel_include_draw(j);
          a_randu = arma::log(arma::randu<arma::vec>(1));
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) >= a_bvs_lprior_1(a_varsel_pos)){continue;}
          if (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) >= a_bvs_lprior_0(a_varsel_pos)){continue;}
          if ((a_lambda(a_varsel_pos, a_varsel_pos) == 1 && a_randu(0) < a_bvs_lprior_1(a_varsel_pos)) || (a_lambda(a_varsel_pos, a_varsel_pos) == 0 && a_randu(0) < a_bvs_lprior_0(a_varsel_pos))){
            a_theta0 = a_AG;
            a_theta1 = a_AG;
            a_theta0.row(a_varsel_pos) = 0;
            a_theta1.row(a_varsel_pos) = gamma.row(a_varsel_pos);
            a_bvs_l0_res = yvec - z * a_theta0;
            a_bvs_l1_res = yvec - z * a_theta1;
            a_l0 = -arma::as_scalar(trans(a_bvs_l0_res) * diag_sigma_i * a_bvs_l0_res) / 2 + arma::as_scalar(a_bvs_lprior_0(a_varsel_pos));
            a_l1 = -arma::as_scalar(trans(a_bvs_l1_res) * diag_sigma_i * a_bvs_l1_res) / 2 + arma::as_scalar(a_bvs_lprior_1(a_varsel_pos));
            a_bayes = a_l1 - a_l0;
            a_bayes_rand = arma::as_scalar(arma::log(arma::randu<arma::vec>(1)));
            if (a_bayes >= a_bayes_rand){
              a_lambda(a_varsel_pos, a_varsel_pos) = 1;
            } else {
              a_lambda(a_varsel_pos, a_varsel_pos) = 0;
            }
          }
        }
        gamma = a_lambda * gamma;
      }
    }
    
    if (use_rr) {
      // Reparameterise alpha
      alpha = arma::reshape(gamma.rows(0, n_alpha - 1), k, r);
      Alpha = alpha * arma::inv(arma::sqrtmat_sympd(alpha.t() * alpha));
      
      // Draw Beta
      y_beta = yvec - z.cols(n_alpha, n_tot - 1) * gamma.rows(n_alpha, n_tot - 1);
      for (int i = 0; i < tt; i++){
        z_beta.rows(i * k, (i + 1) * k - 1) = arma::kron(Alpha, arma::trans(w.col(i)));
      }
      
      post_beta_Vi = arma::solve(arma::kron(Alpha.t() * sigma_i * Alpha, coint_v_i * p_tau_i) + arma::trans(z_beta) * diag_sigma_i * z_beta, diag_beta);
      post_beta_mu = post_beta_Vi * (arma::trans(z_beta) * diag_sigma_i * y_beta);
      Beta = arma::reshape(arma::mvnrnd(post_beta_mu, post_beta_Vi), n_w, r);
      
      // Final cointegration values
      BB_sqrt = arma::sqrtmat_sympd(arma::trans(Beta) * Beta);
      alpha = Alpha * BB_sqrt;
      beta = Beta * arma::inv(BB_sqrt);
      pi = alpha * beta.t();
    } 
    
    // Store draws
    if (draw >= burnin) {
      pos_draw = draw - burnin;
      
      draws_a.col(pos_draw) = gamma;
      if (use_rr) {
        draws_alpha.col(pos_draw) = arma::vectorise(alpha);
        draws_beta.col(pos_draw) = arma::vectorise(beta);
        draws_pi.col(pos_draw) = arma::vectorise(pi);
      }
      if (ssvs) {
        draws_a_lambda.col(pos_draw) = a_lambda.col(0);
      }
      if (bvs) {
        draws_a_lambda.col(pos_draw) = a_lambda.diag();
      }
      draws_sigma.col(pos_draw) = vectorise(arma::solve(sigma_i, diag_k));
      if (vs_covar) {
        if (ssvs) {
          draws_cov_lambda.col(pos_draw) = psi_lambda;
        }
        if (bvs) {
          draws_cov_lambda.col(pos_draw) = cov_lambda.diag();
        }
      }
    }
  } // End loop
  
  // Rcpp::List bla = Rcpp::List::create(Rcpp::Named("test") = sigma_i);
  // return bla;
  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a0") = R_NilValue,
                                             Rcpp::Named("alpha") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
                                             Rcpp::Named("pi") = R_NilValue,
                                             Rcpp::Named("pi_exo") = R_NilValue,
                                             Rcpp::Named("pi_d") = R_NilValue,
                                             Rcpp::Named("gamma") = R_NilValue,
                                             Rcpp::Named("upsilon") = R_NilValue,
                                             Rcpp::Named("c") = R_NilValue,
                                             Rcpp::Named("sigma") = R_NilValue);
  
  if (use_rr) {
    posteriors["alpha"] = Rcpp::wrap(draws_alpha);
    posteriors["beta"] = Rcpp::wrap(draws_beta);
    posteriors["pi"] = Rcpp::wrap(draws_pi.rows(0, k * k - 1));
    if (m > 0) {
      posteriors["pi_exo"] = Rcpp::wrap(draws_pi.rows(k * k, k * (k + m) - 1));
    }
    if (n_r > 0) {
      posteriors["pi_d"] = Rcpp::wrap(draws_pi.rows(k * (k + m), k * (k + m + n_r) - 1));
    }
  }
  
  if (n_gamma > 0) {
    if (varsel) {
      posteriors["gamma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a.rows(n_alpha, n_alpha + n_gamma - 1),
                                                          Rcpp::Named("lambda") = draws_a_lambda.rows(n_alpha, n_alpha + n_gamma - 1)));
    } else {
      posteriors["gamma"] = Rcpp::wrap(draws_a.rows(n_alpha, n_alpha + n_gamma - 1));
    }
  }
  if (n_upsilon > 0) {
    if (varsel) {
      posteriors["upsilon"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a.rows(n_alpha + n_gamma, n_alpha + n_gamma + n_upsilon - 1),
                                                            Rcpp::Named("lambda") = draws_a_lambda.rows(n_alpha + n_gamma, n_alpha + n_gamma + n_upsilon - 1)));
    } else {
      posteriors["upsilon"] = Rcpp::wrap(draws_a.rows(n_alpha + n_gamma, n_alpha + n_gamma + n_upsilon - 1));
    }
  }
  if (n_ur > 0) {
    if (varsel) {
      posteriors["c"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a.rows(n_alpha + n_gamma + n_upsilon, n_alpha + n_gamma + n_upsilon + n_c_ur - 1),
                                                      Rcpp::Named("lambda") = draws_a_lambda.rows(n_alpha + n_gamma + n_upsilon, n_alpha + n_gamma + n_upsilon + n_c_ur - 1)));
    } else {
      posteriors["c"] = Rcpp::wrap(draws_a.rows(n_alpha + n_gamma + n_upsilon, n_alpha + n_gamma + n_upsilon + n_c_ur - 1));
    }
  }
  if (n_a0 > 0) {
    if (varsel) {
      posteriors["a0"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_a.rows(n_alpha + n_gamma + n_upsilon + n_c_ur, n_alpha + n_gamma + n_upsilon + n_c_ur + n_a0 - 1),
                                                       Rcpp::Named("lambda") = draws_a_lambda.rows(n_alpha + n_gamma + n_upsilon + n_c_ur, n_alpha + n_gamma + n_upsilon + n_c_ur + n_a0 - 1)));
    } else {
      posteriors["a0"] = Rcpp::wrap(draws_a.rows(n_alpha + n_gamma + n_upsilon + n_c_ur, n_alpha + n_gamma + n_upsilon + n_c_ur + n_a0 - 1));
    }
  }
  
  if (vs_covar) {
    posteriors["sigma"] = Rcpp::wrap(Rcpp::List::create(Rcpp::Named("coeffs") = draws_sigma,
                                                        Rcpp::Named("lambda") = draws_cov_lambda));
  } else {
    posteriors["sigma"] = Rcpp::wrap(draws_sigma);
  }
  
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("data") = object["data"],
                                         Rcpp::Named("model") = object["model"],
                                         Rcpp::Named("priors") = object["priors"],
                                         Rcpp::Named("posteriors") = posteriors);
  
  return result;
}