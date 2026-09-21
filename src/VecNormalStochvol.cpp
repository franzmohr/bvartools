#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_normal_stochvol.h"

// VEC with a normal prior on the non-cointegration coefficients, the
// cointegration space prior on beta and stochastic volatility in the errors,
// optionally with a constant covariance block. The numerics are the vendored
// BayesTS core; what is here is the translation between the R model object and
// the core's structs, plus the entry points add_posterior_forecasts() and
// add_posterior_loglik() dispatch to. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecNormalStochvolInput read_input(const Rcpp::List &object) {

  bayests::VecNormalStochvolInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model, "sv+covar");

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
      read_forecast_regressors(forecast, input.spec.k, input.forecast.x);
    }
    // What the horizon realised, where the model carries it: one row per
    // period and one column per variable, in levels, which is what
    // /data/test/y of a model file holds and what add_forecast_errors() leaves
    // behind. Absent from every model that is forecast rather than scored, and
    // left for the core to complain about only if a score is actually asked
    // for.
    if (has(data, "test")) {
      const Rcpp::List test = data["test"];
      read_mat_if_present(test, "y", input.test.y);
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
      input.varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }
  if (has(priors, "beta")) {
    input.beta_prior = read_coint_space_prior_constant(priors["beta"]);
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
  read_vec_if_present(initial, "beta", input.initial.beta);
  read_vec_if_present(initial, "psi", input.initial.psi);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);
  read_mat_if_present(initial, "h", input.initial.h);
  read_vec_if_present(initial, "h_init", input.initial.h_init);

  return input;
}

/// The coefficients are constant here and the precision is not, so everything
/// but `u_sigma_inv` is read the same way whichever entry point asked. `beta` is
/// among them: without it `a` carries only the loadings, and both the forecast
/// and the log likelihood rebuild the loadings' regressors from the
/// cointegration matrix.
bayests::VecNormalStochvolDraws read_coefficient_draws(const Rcpp::List &posterior) {

  bayests::VecNormalStochvolDraws draws;

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

  return draws;
}

/// The volatility starts from its last in-sample value, so the forecast wants
/// that period of the precision path alone -- and, when it simulates the
/// volatility forward, the same period of its diagonal and the variance of the
/// log-volatility innovations to step it by.
bayests::VecNormalStochvolDraws read_draws_for_forecast(
    const Rcpp::List &object, const bayests::VecNormalStochvolInput &input) {

  if (!has(object, "posterior")) {
    return bayests::VecNormalStochvolDraws();
  }

  const Rcpp::List posterior = object["posterior"];
  bayests::VecNormalStochvolDraws draws = read_coefficient_draws(posterior);
  const arma::uword k = static_cast<arma::uword>(input.spec.k);

  if (has(posterior, "u_sigma_inv")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs",
                                      input.train.periods(input.spec.k), k * k,
                                      draws.u_sigma_inv);
  }

  if (input.spec.forecast_states == bayests::ForecastStates::simulate) {
    if (has(posterior, "u_omega_inv")) {
      read_draws_last_period_if_present(Rcpp::List(posterior["u_omega_inv"]), "coeffs",
                                        input.train.periods(input.spec.k), k, draws.u_omega_inv);
    }
    if (has(posterior, "u_sigma_inv")) {
      read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "sigma", draws.h_sigma);
    }
  }

  return draws;
}

/// Every period is scored under its own precision, so the log likelihood wants
/// the whole path.
bayests::VecNormalStochvolDraws read_draws_for_loglik(const Rcpp::List &object) {

  if (!has(object, "posterior")) {
    return bayests::VecNormalStochvolDraws();
  }

  const Rcpp::List posterior = object["posterior"];
  bayests::VecNormalStochvolDraws draws = read_coefficient_draws(posterior);

  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// A member left NULL is how "the model did not have that part" is expressed on
/// the R side, which is the convention the core expresses with an empty matrix.
Rcpp::List write_draws(const bayests::VecNormalStochvolDraws &draws) {

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
  // `sigma` is the variance of the log-volatility innovations, as the VAR
  // models keep it.
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv),
                                                 Rcpp::Named("sigma") = draws_to_r(draws.h_sigma));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VecNormalStochvolCoefficients)]]
Rcpp::List VecNormalStochvolCoefficients(Rcpp::List object) {

  const bayests::VecNormalStochvolInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecNormalStochvolDraws draws =
    bayests::VecNormalStochvolSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws),
                            Rcpp::Named("warnings") = reporter.warnings());
}

// [[Rcpp::export(.VecNormalStochvolForecasts)]]
Rcpp::List VecNormalStochvolForecasts(Rcpp::List object) {

  const bayests::VecNormalStochvolInput input = read_input(object);
  const bayests::VecNormalStochvolDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecNormalStochvolSampler().forecast(input, draws, reporter);

  return with_forecast_member(object, "forecasts", Rcpp::wrap(draws_to_r(forecast.values)));
}

// [[Rcpp::export(.VecNormalStochvolLogLik)]]
Rcpp::List VecNormalStochvolLogLik(Rcpp::List object) {

  const bayests::VecNormalStochvolInput input = read_input(object);
  const bayests::VecNormalStochvolDraws draws = read_draws_for_loglik(object);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecNormalStochvolSampler().log_likelihood(input, draws);

  return with_posterior_element(object, "loglik", Rcpp::wrap(loglik));
}

// [[Rcpp::export(.VecNormalStochvolScore)]]
Rcpp::List VecNormalStochvolScore(Rcpp::List object) {

  const bayests::VecNormalStochvolInput input = read_input(object);
  // The same slice a forecast reads: what moves with time is taken at the last
  // in-sample period, which is where carrying it forward over the scored
  // periods starts, together with the state variances each step is drawn with.
  const bayests::VecNormalStochvolDraws draws = read_draws_for_forecast(object, input);

  // Draws by scored periods already, the same orientation the pointwise log
  // likelihood comes back in, so there is no transpose at this boundary either.
  const arma::mat score = bayests::VecNormalStochvolSampler().predictive_log_density(input, draws);

  return with_forecast_member(object, "loglik", Rcpp::wrap(score));
}

/*** R

data("us_macrodata")

object <- create_bvecmodel(data = us_macrodata,
                           p = 2, const = "unrestricted",
                           r = 1,
                           error = "sv+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = .01,
                                  state_variance = .05, offset = .0001))

object <- add_initial_values(object)

object <- .VecNormalStochvolCoefficients(object)
object <- .VecNormalStochvolLogLik(object)

object <- add_forecast_input(object, n_ahead = 10)
object <- .VecNormalStochvolForecasts(object)

*/
