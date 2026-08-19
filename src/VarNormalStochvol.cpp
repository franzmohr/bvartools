#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_normal_stochvol.h"

// VAR with a normal prior on the coefficients and stochastic volatility in the
// errors, optionally with a covariance block. The numerics are the vendored
// BayesTS core; what is here is the translation between the R model object and
// the core's structs. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VarNormalStochvolInput read_input(const Rcpp::List &object) {

  bayests::VarNormalStochvolInput input;

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

  const Rcpp::List initial = has(object, "initial") ? Rcpp::List(object["initial"]) : Rcpp::List();
  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();

  // BVS is the only scheme this model implements. An SSVS object is left unread
  // here and rejected by validate(), which can say why -- reading it would fail
  // first, on mixture components the object does not carry.
  const bool bvs = input.spec.varsel == bayests::VarSelection::bvs;

  if (has(priors, "a")) {
    const Rcpp::List prior_a = priors["a"];
    input.a_prior = read_normal_prior(prior_a);
    if (bvs) {
      input.a_varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }
  if (has(priors, "psi")) {
    const Rcpp::List prior_psi = priors["psi"];
    input.psi_prior = read_normal_prior(prior_psi);
    if (bvs) {
      input.psi_varsel_prior = read_varsel_prior(prior_psi, input.spec.varsel);
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

  read_vec_if_present(initial, "a", input.initial.a);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);
  read_vec_if_present(initial, "psi", input.initial.psi);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);
  read_mat_if_present(initial, "h", input.initial.h);
  read_vec_if_present(initial, "h_init", input.initial.h_init);

  return input;
}

bayests::VarNormalStochvolDraws read_coefficient_draws(const Rcpp::List &posterior) {

  bayests::VarNormalStochvolDraws draws;

  if (has(posterior, "a")) {
    const Rcpp::List posterior_a = posterior["a"];
    read_draws_if_present(posterior_a, "coeffs", draws.a);
    read_draws_if_present(posterior_a, "lambda", draws.a_lambda);
  }
  if (has(posterior, "psi")) {
    const Rcpp::List posterior_psi = posterior["psi"];
    read_draws_if_present(posterior_psi, "coeffs", draws.psi);
    read_draws_if_present(posterior_psi, "lambda", draws.psi_lambda);
  }

  return draws;
}

/// The volatility is held at its last in-sample value, so the forecast wants
/// that period of the precision path alone.
bayests::VarNormalStochvolDraws read_draws_for_forecast(const Rcpp::List &object,
                                                        const bayests::VarNormalStochvolInput &input) {

  if (!has(object, "posterior")) {
    return bayests::VarNormalStochvolDraws();
  }

  const Rcpp::List posterior = object["posterior"];
  bayests::VarNormalStochvolDraws draws = read_coefficient_draws(posterior);

  if (has(posterior, "u_sigma_inv")) {
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs",
                                      input.train.periods(input.spec.k), k * k,
                                      draws.u_sigma_inv);
  }

  return draws;
}

/// Every period is scored under its own precision, so the log likelihood wants
/// the whole path.
bayests::VarNormalStochvolDraws read_draws_for_loglik(const Rcpp::List &object) {

  if (!has(object, "posterior")) {
    return bayests::VarNormalStochvolDraws();
  }

  const Rcpp::List posterior = object["posterior"];
  bayests::VarNormalStochvolDraws draws = read_coefficient_draws(posterior);

  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

Rcpp::List write_draws(const bayests::VarNormalStochvolDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("psi") = R_NilValue,
                                             Rcpp::Named("u_sigma_inv") = R_NilValue);

  if (draws.has_a()) {
    if (draws.a_lambda.n_elem > 0) {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                           Rcpp::Named("lambda") = draws_to_r(draws.a_lambda));
    } else {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a));
    }
  }

  if (draws.has_psi()) {
    if (draws.psi_lambda.n_elem > 0) {
      posteriors["psi"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.psi),
                                             Rcpp::Named("lambda") = draws_to_r(draws.psi_lambda));
    } else {
      posteriors["psi"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.psi));
    }
  }

  posteriors["u_omega_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_omega_inv));
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VarNormalStochvolCoefficients)]]
Rcpp::List VarNormalStochvolCoefficients(Rcpp::List object) {

  const bayests::VarNormalStochvolInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarNormalStochvolDraws draws =
    bayests::VarNormalStochvolSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VarNormalStochvolForecasts)]]
Rcpp::List VarNormalStochvolForecasts(Rcpp::List object) {

  const bayests::VarNormalStochvolInput input = read_input(object);
  const bayests::VarNormalStochvolDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarNormalStochvolSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VarNormalStochvolLogLik)]]
Rcpp::List VarNormalStochvolLogLik(Rcpp::List object) {

  const bayests::VarNormalStochvolInput input = read_input(object);
  const bayests::VarNormalStochvolDraws draws = read_draws_for_loglik(object);

  const arma::mat loglik = bayests::VarNormalStochvolSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata,
                           p = 0, deterministic = "const",
                           error = "sv+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = .01,
                                  state_variance = .05, offset = .0001))

object <- add_initial_values(object)

object <- .VarNormalStochvolCoefficients(object)
object <- .VarNormalStochvolLogLik(object)

*/
