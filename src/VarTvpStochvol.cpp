#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_tvp_stochvol.h"

// VAR whose coefficients follow a random walk, with stochastic volatility in
// the errors and an optional time-varying covariance block. The numerics are
// the vendored BayesTS core; what is here is the translation between the R
// model object and the core's structs. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VarTvpStochvolInput read_input(const Rcpp::List &object) {

  bayests::VarTvpStochvolInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model, "sv+covar");

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "z", input.train.z);
    }
    if (has(data, "forecast")) {
      const Rcpp::List forecast = data["forecast"];
      read_mat_if_present(forecast, "z", input.forecast.z);
    }
  }

  const arma::uword tt = input.train.y.n_elem > 0 && input.spec.k > 0
                           ? input.train.periods(input.spec.k)
                           : 0;

  const Rcpp::List initial = has(object, "initial") ? Rcpp::List(object["initial"]) : Rcpp::List();
  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();

  // BVS is the only scheme this model implements. An SSVS object is left unread
  // here and rejected by validate(), which can say why.
  const bool bvs = input.spec.varsel == bayests::VarSelection::bvs;

  if (has(priors, "a")) {
    const Rcpp::List prior_a = priors["a"];
    // One R group carries both halves of the state equation: how far the
    // coefficients may drift, and where they start.
    input.a_prior.sigma = read_gamma_prior(prior_a);
    input.a_prior.initial_state = read_normal_prior(prior_a);
    if (bvs) {
      input.a_varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }

  if (has(priors, "psi")) {
    const Rcpp::List prior_psi = priors["psi"];
    input.psi_prior.sigma = read_gamma_prior(prior_psi);
    input.psi_prior.initial_state = read_normal_prior(prior_psi);

    // Selection for the covariance block is declared in its own group, so it
    // can differ from the model's.
    input.psi_varsel = bayests::var_selection_from_string(
      optional_string(prior_psi, "varsel", "none"));

    if (input.psi_varsel == bayests::VarSelection::bvs) {
      input.psi_varsel_prior = read_varsel_prior(prior_psi, input.psi_varsel);
    }
  }

  if (has(priors, "u_sigma")) {
    const Rcpp::List prior_u_sigma = priors["u_sigma"];
    read_vec_if_present(prior_u_sigma, "offset", input.u_sigma_prior.offset);
    read_vec_if_present(prior_u_sigma, "shape", input.u_sigma_prior.state.sigma.shape);
    read_vec_if_present(prior_u_sigma, "rate", input.u_sigma_prior.state.sigma.rate);
    read_vec_if_present(prior_u_sigma, "mu", input.u_sigma_prior.state.initial_state.mu);
    read_mat_if_present(prior_u_sigma, "v_inv", input.u_sigma_prior.state.initial_state.v_inv);
    // A state the sampler redraws every iteration, even though R keeps it next
    // to the prior it is drawn under.
    read_vec_if_present(prior_u_sigma, "sigma", input.initial.h_sigma);
  }

  read_path_if_present(initial, "a", input.train.nparams(), tt, input.initial.a);
  read_mat_if_present(initial, "a_sigma_inv", input.initial.a_sigma_inv);
  read_vec_if_present(initial, "a_init", input.initial.a_init);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);

  read_path_if_present(initial, "psi", static_cast<arma::uword>(input.spec.n_psi()), tt,
                       input.initial.psi);
  read_mat_if_present(initial, "psi_sigma_inv", input.initial.psi_sigma_inv);
  read_vec_if_present(initial, "psi_init", input.initial.psi_init);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);

  read_mat_if_present(initial, "h", input.initial.h);
  read_vec_if_present(initial, "h_init", input.initial.h_init);

  return input;
}

/// Both the coefficients and the precision are held at their last in-sample
/// values, so both are sliced to that period.
bayests::VarTvpStochvolDraws read_draws_for_forecast(const Rcpp::List &object,
                                                     const bayests::VarTvpStochvolInput &input) {

  bayests::VarTvpStochvolDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  const arma::uword tt = input.train.periods(input.spec.k);
  const arma::uword k = static_cast<arma::uword>(input.spec.k);

  if (input.forecast.z.n_cols > 0 && has(posterior, "a")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["a"]), "coeffs", tt,
                                      static_cast<arma::uword>(input.spec.nparams_per_period()),
                                      draws.a);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", tt, k * k,
                                      draws.u_sigma_inv);
  }

  return draws;
}

/// Every period is evaluated under its own coefficients, so `a` is the whole
/// path here; the precision is the one period this model scores everything
/// under.
bayests::VarTvpStochvolDraws read_draws_for_loglik(const Rcpp::List &object,
                                                   const bayests::VarTvpStochvolInput &input) {

  bayests::VarTvpStochvolDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  const arma::uword k = static_cast<arma::uword>(input.spec.k);

  if (has(posterior, "a")) {
    read_draws_if_present(Rcpp::List(posterior["a"]), "coeffs", draws.a);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs",
                                      input.train.periods(input.spec.k), k * k,
                                      draws.u_sigma_inv);
  }

  return draws;
}

Rcpp::List write_draws(const bayests::VarTvpStochvolDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("psi") = R_NilValue,
                                             Rcpp::Named("u_sigma_inv") = R_NilValue);

  if (draws.has_a()) {
    if (draws.a_lambda.n_elem > 0) {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                           Rcpp::Named("sigma") = draws_to_r(draws.a_sigma),
                                           Rcpp::Named("lambda") = draws_to_r(draws.a_lambda));
    } else {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                           Rcpp::Named("sigma") = draws_to_r(draws.a_sigma));
    }
  }

  if (draws.has_psi()) {
    if (draws.psi_lambda.n_elem > 0) {
      posteriors["psi"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.psi),
                                             Rcpp::Named("sigma") = draws_to_r(draws.psi_sigma),
                                             Rcpp::Named("lambda") = draws_to_r(draws.psi_lambda));
    } else {
      posteriors["psi"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.psi),
                                             Rcpp::Named("sigma") = draws_to_r(draws.psi_sigma));
    }
  }

  posteriors["u_omega_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_omega_inv));
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VarTvpStochvolCoefficients)]]
Rcpp::List VarTvpStochvolCoefficients(Rcpp::List object) {

  const bayests::VarTvpStochvolInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarTvpStochvolDraws draws =
    bayests::VarTvpStochvolSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VarTvpStochvolForecasts)]]
Rcpp::List VarTvpStochvolForecasts(Rcpp::List object) {

  const bayests::VarTvpStochvolInput input = read_input(object);
  const bayests::VarTvpStochvolDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarTvpStochvolSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VarTvpStochvolLogLik)]]
Rcpp::List VarTvpStochvolLogLik(Rcpp::List object) {

  const bayests::VarTvpStochvolInput input = read_input(object);
  const bayests::VarTvpStochvolDraws draws = read_draws_for_loglik(object, input);

  const arma::mat loglik = bayests::VarTvpStochvolSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata,
                           p = 0, deterministic = "const",
                           tvp = TRUE, error = "sv+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1, shape = 3, rate = .0001),
                     sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = .01,
                                  state_variance = .05, offset = .0001))

object <- add_initial_values(object)

object <- .VarTvpStochvolCoefficients(object)
object <- .VarTvpStochvolLogLik(object)

*/
