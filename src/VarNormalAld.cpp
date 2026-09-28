#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_normal_ald.h"

// VAR estimated at a conditional quantile rather than at the conditional mean,
// through the normal scale mixture representation of the asymmetric Laplace
// distribution. The numerics are the vendored BayesTS core; what is here is the
// translation between the R model object and the core's structs. See
// src/core/VENDORED.md.
//
// It takes no covariance block, and at a single quantile it does not forecast:
// the core's validate() refuses both with a reason. A grid of quantiles,
// model$quantiles, is the structural quantile VAR of Chavleishvili and
// Manganelli (2019), which describes the whole conditional distribution and
// forecasts by simulating through it -- what .VarNormalAldForecasts is for.

namespace {

using namespace bayests_r;

bayests::VarNormalAldInput read_input(const Rcpp::List &object) {

  bayests::VarNormalAldInput input;

  const Rcpp::List model = object["model"];
  // No spelling of 'error' means "with a covariance block" for this model, so
  // spec.covar stays false and an object that asked for one is refused by
  // validate() rather than silently estimated without it.
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "z", input.train.z);
    }
    // Read for a grid of quantiles, and refused by validate() for a single
    // quantile together with the non-zero horizon add_forecast_input() sets.
    if (has(data, "forecast")) {
      const Rcpp::List forecast = data["forecast"];
      read_forecast_regressors(forecast, input.spec.k, input.forecast.x);
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
  if (has(priors, "u_scale")) {
    input.u_scale_prior = read_gamma_prior(Rcpp::List(priors["u_scale"]));
  }

  read_vec_if_present(initial, "a", input.initial.a);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);
  // The latent scales are tt x k, one column per equation, and the sampler
  // redraws every element of them each sweep.
  read_mat_if_present(initial, "w", input.initial.w);
  read_vec_if_present(initial, "u_scale", input.initial.u_scale);

  read_extensions(object, input);
  return input;
}

/// The log likelihood and the forecast of a grid want the coefficients and the
/// scale and nothing else: the asymmetric Laplace density is closed form and
/// marginal of the latent scales, and a grid simulates from the distribution
/// its quantiles describe. The precision path is read for its column count
/// alone -- iterations() counts it -- and one period is enough for that, which
/// saves carrying k * k * tt numbers per draw across the boundary.
bayests::VarNormalAldDraws read_draws(const Rcpp::List &object,
                                                 const bayests::VarNormalAldInput &input) {

  bayests::VarNormalAldDraws draws;

  if (!has(object, "posterior")) {
    read_draw_extensions(object, draws);
  return draws;
  }

  const Rcpp::List posterior = object["posterior"];

  if (has(posterior, "a")) {
    const Rcpp::List posterior_a = posterior["a"];
    read_draws_if_present(posterior_a, "coeffs", draws.a);
    read_draws_if_present(posterior_a, "lambda", draws.a_lambda);
  }
  if (has(posterior, "u_scale")) {
    read_draws_if_present(Rcpp::List(posterior["u_scale"]), "coeffs", draws.u_scale);
  }
  if (has(posterior, "u_sigma_inv")) {
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    read_draws_last_period_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs",
                                      input.train.periods(input.spec.k), k * k,
                                      draws.u_sigma_inv);
  }

  read_draw_extensions(object, draws);
  return draws;
}

Rcpp::List write_draws(const bayests::VarNormalAldDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("u_sigma_inv") = R_NilValue);

  if (draws.has_a()) {
    if (draws.a_lambda.n_elem > 0) {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a),
                                           Rcpp::Named("lambda") = draws_to_r(draws.a_lambda));
    } else {
      posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a));
    }
  }

  // The scale of the asymmetric Laplace, one per equation. The latent scales it
  // multiplies are not returned: they are k * tt numbers of pure nuisance per
  // draw, and u_omega_inv together with u_scale recovers them.
  posteriors["u_scale"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_scale));
  // A grid of quantiles keeps no latent scales -- each level's are a nuisance
  // of its own chain -- so the two are left NULL and dropped on the R side.
  if (draws.u_omega_inv.n_elem > 0) {
    posteriors["u_omega_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_omega_inv));
  }
  if (draws.u_sigma_inv.n_elem > 0) {
    posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));
  }

  write_draw_extensions(posteriors, draws);
  return posteriors;
}

} // namespace

// [[Rcpp::export(.VarNormalAldCoefficients)]]
Rcpp::List VarNormalAldCoefficients(Rcpp::List object) {

  const bayests::VarNormalAldInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarNormalAldDraws draws =
    bayests::VarNormalAldSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws),
                            Rcpp::Named("warnings") = reporter.warnings());
}

// [[Rcpp::export(.VarNormalAldForecasts)]]
Rcpp::List VarNormalAldForecasts(Rcpp::List object) {

  const bayests::VarNormalAldInput input = read_input(object);
  const bayests::VarNormalAldDraws draws = read_draws(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarNormalAldSampler().forecast(input, draws, reporter);

  return with_forecast_member(object, "forecasts", Rcpp::wrap(draws_to_r(forecast.values)));
}

// [[Rcpp::export(.VarNormalAldLogLik)]]
Rcpp::List VarNormalAldLogLik(Rcpp::List object) {

  const bayests::VarNormalAldInput input = read_input(object);
  const bayests::VarNormalAldDraws draws = read_draws(object, input);

  const arma::mat loglik = bayests::VarNormalAldSampler().log_likelihood(input, draws);

  return with_posterior_element(object, "loglik", Rcpp::wrap(loglik));
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata,
                           p = 1, deterministic = "const",
                           error = "ald", quantile = 0.25,
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     sigma = list(shape = 3, rate = .01))

object <- add_initial_values(object)

object <- .VarNormalAldCoefficients(object)
object <- .VarNormalAldLogLik(object)

*/
