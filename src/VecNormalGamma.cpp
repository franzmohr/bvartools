#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_normal_gamma.h"

// VEC with a normal prior on the non-cointegration coefficients, the
// cointegration space prior on beta and independent gamma priors on the error
// precisions, optionally with a constant covariance block. The numerics are the
// vendored BayesTS core; what is here is the translation between the R model
// object and the core's structs, plus the entry points
// add_posterior_forecasts() and add_posterior_loglik() dispatch to. See
// src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecNormalGammaInput read_input(const Rcpp::List &object) {

  bayests::VecNormalGammaInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model, "gamma+covar");

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      // The cointegration regressors, which a VAR has no counterpart for.
      read_mat_if_present(train, "w", input.train.w);
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
      // R keeps the selection prior alongside the normal one rather than in a
      // group of its own.
      input.varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }
  if (has(priors, "beta")) {
    input.beta_prior = read_coint_space_prior_constant(priors["beta"]);
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
  read_vec_if_present(initial, "beta", input.initial.beta);
  read_vec_if_present(initial, "psi", input.initial.psi);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);
  // R calls the starting precision u_omega_inv; the core calls the same matrix
  // u_sigma_inv, and the sampler opens by copying one into the other.
  read_mat_if_present(initial, "u_omega_inv", input.initial.u_sigma_inv);

  return input;
}

/// Nothing in this model moves with time, so the forecast and the log
/// likelihood read every draw as it stands. `beta` is read for both: without it
/// `a` carries only the loadings, and each rebuilds the loadings' regressors
/// from the cointegration matrix.
bayests::VecNormalGammaDraws read_draws(const Rcpp::List &object) {

  bayests::VecNormalGammaDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  if (has(posterior, "a")) {
    const Rcpp::List posterior_a = posterior["a"];
    read_draws_if_present(posterior_a, "coeffs", draws.a);
    read_draws_if_present(posterior_a, "lambda", draws.a_lambda);
  }
  if (has(posterior, "beta")) {
    read_draws_if_present(Rcpp::List(posterior["beta"]), "coeffs", draws.beta);
  }
  if (has(posterior, "psi")) {
    const Rcpp::List posterior_psi = posterior["psi"];
    read_draws_if_present(posterior_psi, "coeffs", draws.psi);
    read_draws_if_present(posterior_psi, "lambda", draws.psi_lambda);
  }
  if (has(posterior, "u_omega_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_omega_inv"]), "coeffs", draws.u_omega_inv);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// A member left NULL is how "the model did not have that part" is expressed on
/// the R side, which is the convention the core expresses with an empty matrix.
Rcpp::List write_draws(const bayests::VecNormalGammaDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
                                             Rcpp::Named("psi") = R_NilValue,
                                             Rcpp::Named("u_sigma_inv") = R_NilValue);

  if (draws.has_a()) {
    if (draws.has_lambda()) {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                           Rcpp::Named("lambda") = draws_to_r(draws.a_lambda));
    } else {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a));
    }
  }

  if (draws.has_beta()) {
    posteriors["beta"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.beta));
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

// [[Rcpp::export(.VecNormalGammaCoefficients)]]
Rcpp::List VecNormalGammaCoefficients(Rcpp::List object) {

  const bayests::VecNormalGammaInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecNormalGammaDraws draws =
    bayests::VecNormalGammaSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VecNormalGammaForecasts)]]
Rcpp::List VecNormalGammaForecasts(Rcpp::List object) {

  const bayests::VecNormalGammaInput input = read_input(object);
  const bayests::VecNormalGammaDraws draws = read_draws(object);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecNormalGammaSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VecNormalGammaLogLik)]]
Rcpp::List VecNormalGammaLogLik(Rcpp::List object) {

  const bayests::VecNormalGammaInput input = read_input(object);
  const bayests::VecNormalGammaDraws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecNormalGammaSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvecmodel(data = us_macrodata,
                           p = 2, const = "unrestricted",
                           r = 1,
                           error = "gamma+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(shape = 3, rate = 1))

object <- add_initial_values(object)

object <- .VecNormalGammaCoefficients(object)
object <- .VecNormalGammaLogLik(object)

object <- add_forecast_input_data(object, n_ahead = 10)
object <- .VecNormalGammaForecasts(object)

*/
