#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_tvp_stochvol.h"

// VEC whose loadings, lagged differences and cointegration vectors all follow a
// random walk, with stochastic volatility in the errors and an optional
// time-varying covariance block. The numerics are the vendored BayesTS core;
// what is here is the translation between the R model object and the core's
// structs, plus the entry points add_posterior_forecasts() and
// add_posterior_loglik() dispatch to. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecTvpStochvolInput read_input(const Rcpp::List &object) {

  bayests::VecTvpStochvolInput input;

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

  // Every path is stored flat in R and wanted one column per period, so the
  // number of periods has to be known before any of them can be read.
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
    // The non-centred parameterisation, in place of shape and rate; validate()
    // refuses a block that carries both.
    read_vec_if_present(prior_a, "omega_v", input.a_prior.omega_v);
    if (bvs) {
      input.a_varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }

  if (has(priors, "beta")) {
    input.beta_prior = read_coint_space_prior_tvp(priors["beta"]);
  }

  if (has(priors, "psi")) {
    const Rcpp::List prior_psi = priors["psi"];
    input.psi_prior.sigma = read_gamma_prior(prior_psi);
    input.psi_prior.initial_state = read_normal_prior(prior_psi);
    read_vec_if_present(prior_psi, "omega_v", input.psi_prior.omega_v);

    // Selection for the covariance block is declared in its own group, so it can
    // differ from the model's.
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
    read_vec_if_present(prior_u_sigma, "omega_v", input.u_sigma_prior.state.omega_v);
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

  read_path_if_present(initial, "beta", static_cast<arma::uword>(input.spec.n_beta()), tt,
                       input.initial.beta);
  read_vec_if_present(initial, "beta_init", input.initial.beta_init);

  read_path_if_present(initial, "psi", static_cast<arma::uword>(input.spec.n_psi()), tt,
                       input.initial.psi);
  read_mat_if_present(initial, "psi_sigma_inv", input.initial.psi_sigma_inv);
  read_vec_if_present(initial, "psi_init", input.initial.psi_init);
  read_vec_if_present(initial, "psi_lambda", input.initial.psi_lambda);

  read_mat_if_present(initial, "h", input.initial.h);
  read_vec_if_present(initial, "h_init", input.initial.h_init);

  return input;
}

/// Every period is evaluated under its own coefficients and its own precision,
/// so both paths and the precision are read whole.
bayests::VecTvpStochvolDraws read_draws_for_loglik(const Rcpp::List &object,
                                                   const bayests::VecTvpStochvolInput &input) {

  bayests::VecTvpStochvolDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];

  if (has(posterior, "a")) {
    read_draws_if_present(Rcpp::List(posterior["a"]), "coeffs", draws.a);
  }
  if (has(posterior, "beta")) {
    read_draws_if_present(Rcpp::List(posterior["beta"]), "coeffs", draws.beta);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// The forecast starts from where the sample ended, so the coefficients, the
/// cointegration vectors and the precision are all sliced to that period. The
/// two paths' widths differ and neither is `z.n_cols`: `a` counts off the spec,
/// since a VEC stores the loadings in it, and `beta` off the cointegration
/// dimensions.
bayests::VecTvpStochvolDraws read_draws_for_forecast(const Rcpp::List &object,
                                                      const bayests::VecTvpStochvolInput &input) {

  bayests::VecTvpStochvolDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  const arma::uword tt = input.train.periods(input.spec.k);
  const arma::uword k = static_cast<arma::uword>(input.spec.k);
  const arma::uword n_a = static_cast<arma::uword>(input.spec.nparams_per_period_vec());
  const arma::uword n_beta = static_cast<arma::uword>(input.spec.n_beta());

  if (n_a > 0 && has(posterior, "a")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["a"]), "coeffs", tt, n_a, draws.a);
  }
  if (n_beta > 0 && has(posterior, "beta")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["beta"]), "coeffs", tt, n_beta,
                                      draws.beta);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", tt, k * k,
                                      draws.u_sigma_inv);
  }

  // What simulating the states forward reads on top of the period they start
  // from: how far each random walk moves per period, which coefficients
  // selection left out, rho where the chain drew it, and the last period of Psi
  // and of the log-volatilities the precision is rebuilt from. A held forecast
  // reads none of it.
  if (input.spec.forecast_states == bayests::ForecastStates::simulate) {
    if (n_a > 0 && has(posterior, "a")) {
      const Rcpp::List block = posterior["a"];
      read_draws_if_present(block, "sigma", draws.a_sigma);
      read_draws_if_present(block, "lambda", draws.a_lambda);
    }
    if (n_beta > 0 && has(posterior, "beta")) {
      read_draws_if_present(Rcpp::List(posterior["beta"]), "rho", draws.rho);
    }
    if (input.use_psi() && has(posterior, "psi")) {
      const Rcpp::List block = posterior["psi"];
      read_draws_last_period_if_present(block, "coeffs", tt, k * k, draws.psi);
      read_draws_if_present(block, "sigma", draws.psi_sigma);
      read_draws_if_present(block, "lambda", draws.psi_lambda);
    }
    if (has(posterior, "u_omega_inv")) {
      read_draws_last_period_if_present(Rcpp::List(posterior["u_omega_inv"]), "coeffs", tt, k,
                                        draws.u_omega_inv);
    }
    if (has(posterior, "u_sigma_inv")) {
      read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "sigma", draws.h_sigma);
    }
  }

  return draws;
}

/// A member left NULL is how "the model did not have that part" is expressed on
/// the R side, which is the convention the core expresses with an empty matrix.
Rcpp::List write_draws(const bayests::VecTvpStochvolDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
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
    posteriors["a"] = with_noncentred(posteriors["a"], draws.a_noncentred);
  }

  // The cointegration path. Without it `a` carries only the loadings, so nothing
  // downstream could reconstruct Pi -- and neither the forecast nor the log
  // likelihood, both of which rebuild the loadings' regressors from it, could be
  // computed at all.
  if (draws.has_beta()) {
    // `rho` only where the model file put a prior on it. Held fixed it is a
    // hyperparameter the caller already has in priors$beta$rho, and handing a
    // column of the same number back would read as a posterior.
    if (draws.has_rho()) {
      posteriors["beta"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.beta),
                                              Rcpp::Named("rho") = draws_to_r(draws.rho));
    } else {
      posteriors["beta"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.beta));
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
    posteriors["psi"] = with_noncentred(posteriors["psi"], draws.psi_noncentred);
  }

  posteriors["u_omega_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_omega_inv));
  // `sigma` is the variance of the log-volatility innovations, as the VAR
  // models keep it.
  posteriors["u_sigma_inv"] = with_noncentred(
    Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv),
                       Rcpp::Named("sigma") = draws_to_r(draws.h_sigma)),
    draws.h_noncentred);

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VecTvpStochvolCoefficients)]]
Rcpp::List VecTvpStochvolCoefficients(Rcpp::List object) {

  const bayests::VecTvpStochvolInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecTvpStochvolDraws draws =
    bayests::VecTvpStochvolSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws),
                            Rcpp::Named("warnings") = reporter.warnings());
}

// [[Rcpp::export(.VecTvpStochvolForecasts)]]
Rcpp::List VecTvpStochvolForecasts(Rcpp::List object) {

  const bayests::VecTvpStochvolInput input = read_input(object);
  const bayests::VecTvpStochvolDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecTvpStochvolSampler().forecast(input, draws, reporter);

  return with_forecast_member(object, "forecasts", Rcpp::wrap(draws_to_r(forecast.values)));
}

// [[Rcpp::export(.VecTvpStochvolLogLik)]]
Rcpp::List VecTvpStochvolLogLik(Rcpp::List object) {

  const bayests::VecTvpStochvolInput input = read_input(object);
  const bayests::VecTvpStochvolDraws draws = read_draws_for_loglik(object, input);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecTvpStochvolSampler().log_likelihood(input, draws);

  return with_posterior_element(object, "loglik", Rcpp::wrap(loglik));
}

// [[Rcpp::export(.VecTvpStochvolScore)]]
Rcpp::List VecTvpStochvolScore(Rcpp::List object) {

  const bayests::VecTvpStochvolInput input = read_input(object);
  // The same slice a forecast reads: what moves with time is taken at the last
  // in-sample period, which is where carrying it forward over the scored
  // periods starts, together with the state variances each step is drawn with.
  const bayests::VecTvpStochvolDraws draws = read_draws_for_forecast(object, input);

  // Draws by scored periods already, the same orientation the pointwise log
  // likelihood comes back in, so there is no transpose at this boundary either.
  const arma::mat score = bayests::VecTvpStochvolSampler().predictive_log_density(input, draws);

  return with_forecast_member(object, "loglik", Rcpp::wrap(score));
}

/*** R

data("us_macrodata")

object <- create_bvecmodel(data = us_macrodata,
                           p = 2, const = "unrestricted",
                           r = 1, tvp = TRUE,
                           error = "sv+covar",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1, shape = 3, rate = .0001),
                     coint = list(rho = .999),
                     sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = .01,
                                  state_variance = .05, offset = .0001))

object <- add_initial_values(object)

object <- .VecTvpStochvolCoefficients(object)
object <- .VecTvpStochvolLogLik(object)

object <- add_forecast_input(object, n_ahead = 10)
object <- .VecTvpStochvolForecasts(object)

*/
