#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_normal_gamma.h"

// VAR with a normal prior on the coefficients and independent gamma priors on
// the error precisions, optionally with a constant covariance block. The
// numerics are the vendored BayesTS core; what is here is the translation
// between the R model object and the core's structs. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VarNormalGammaInput read_input(const Rcpp::List &object) {

  bayests::VarNormalGammaInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model, "gamma+covar");

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

  if (has(priors, "a")) {
    const Rcpp::List prior_a = priors["a"];
    input.a_prior = read_normal_prior(prior_a);
    if (input.spec.uses_varsel()) {
      input.a_varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }
  if (has(priors, "psi")) {
    const Rcpp::List prior_psi = priors["psi"];
    input.psi_prior = read_normal_prior(prior_psi);
    if (input.spec.uses_varsel()) {
      input.psi_varsel_prior = read_varsel_prior(prior_psi, input.spec.varsel);
    }
  }
  if (has(priors, "u_sigma")) {
    input.u_sigma_prior = read_gamma_prior(priors["u_sigma"]);
  }

  read_vec_if_present(initial, "a", input.initial.a);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);
  read_vec_if_present(initial, "psi", input.initial.psi);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);
  // R calls the starting precision u_omega_inv; the core calls the same matrix
  // u_sigma_inv, and both samplers open by copying one into the other.
  read_mat_if_present(initial, "u_omega_inv", input.initial.u_sigma_inv);

  return input;
}

/// Neither the coefficients nor the precision move with time here, so both the
/// forecast and the log likelihood read every draw as it stands. The structural
/// split of `a` is the sampler's business, not this layer's.
bayests::VarNormalGammaDraws read_draws(const Rcpp::List &object) {

  bayests::VarNormalGammaDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
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
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }
  if (has(posterior, "u_omega_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_omega_inv"]), "coeffs", draws.u_omega_inv);
  }

  return draws;
}

Rcpp::List write_draws(const bayests::VarNormalGammaDraws &draws) {

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

// [[Rcpp::export(.VarNormalGammaCoefficients)]]
Rcpp::List VarNormalGammaCoefficients(Rcpp::List object) {

  const bayests::VarNormalGammaInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarNormalGammaDraws draws =
    bayests::VarNormalGammaSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VarNormalGammaForecasts)]]
Rcpp::List VarNormalGammaForecasts(Rcpp::List object) {

  const bayests::VarNormalGammaInput input = read_input(object);
  const bayests::VarNormalGammaDraws draws = read_draws(object);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarNormalGammaSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VarNormalGammaLogLik)]]
Rcpp::List VarNormalGammaLogLik(Rcpp::List object) {

  const bayests::VarNormalGammaInput input = read_input(object);
  const bayests::VarNormalGammaDraws draws = read_draws(object);

  const arma::mat loglik = bayests::VarNormalGammaSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata,
                           p = 0, deterministic = "const",
                           error = "gamma+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     sigma = list(shape = 3, rate = 1))

object <- add_initial_values(object)

object <- .VarNormalGammaCoefficients(object)
object <- .VarNormalGammaLogLik(object)

*/
