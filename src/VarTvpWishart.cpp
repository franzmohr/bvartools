#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_tvp_wishart.h"

// VAR whose coefficients follow a random walk, with a Wishart prior on the
// error precision. The numerics are the vendored BayesTS core; what is here is
// the translation between the R model object and the core's structs. See
// src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VarTvpWishartInput read_input(const Rcpp::List &object) {

  bayests::VarTvpWishartInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

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

  if (has(priors, "u_sigma")) {
    const Rcpp::List prior_u_sigma = priors["u_sigma"];
    input.u_sigma_prior.df = optional_int(prior_u_sigma, "df", 0);
    read_mat_if_present(prior_u_sigma, "scale", input.u_sigma_prior.scale);
  }

  read_path_if_present(initial, "a", input.train.nparams(), tt, input.initial.a);
  read_mat_if_present(initial, "a_sigma_inv", input.initial.a_sigma_inv);
  read_vec_if_present(initial, "a_init", input.initial.a_init);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);
  read_mat_if_present(initial, "u_sigma_inv", input.initial.u_sigma_inv);

  return input;
}

/// The forecast holds the coefficients at their last in-sample value. The
/// precision does not move with time in this model, so it is read whole.
bayests::VarTvpWishartDraws read_draws_for_forecast(const Rcpp::List &object,
                                                    const bayests::VarTvpWishartInput &input) {

  bayests::VarTvpWishartDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];

  // Counted off the forecast regressors, which is what the coefficients drawn
  // per period have to line up with.
  const arma::uword nparams = input.forecast.z.n_cols;
  if (nparams > 0 && has(posterior, "a")) {
    read_draws_last_period_if_present(Rcpp::List(posterior["a"]), "coeffs",
                                      input.train.periods(input.spec.k), nparams, draws.a);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// Every period is evaluated under its own coefficients, so `a` is the whole
/// path here.
bayests::VarTvpWishartDraws read_draws_for_loglik(const Rcpp::List &object) {

  bayests::VarTvpWishartDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  if (has(posterior, "a")) {
    read_draws_if_present(Rcpp::List(posterior["a"]), "coeffs", draws.a);
  }
  if (has(posterior, "u_sigma_inv")) {
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

Rcpp::List write_draws(const bayests::VarTvpWishartDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
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
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VarTvpWishartCoefficients)]]
Rcpp::List VarTvpWishartCoefficients(Rcpp::List object) {

  const bayests::VarTvpWishartInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarTvpWishartDraws draws =
    bayests::VarTvpWishartSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VarTvpWishartForecasts)]]
Rcpp::List VarTvpWishartForecasts(Rcpp::List object) {

  const bayests::VarTvpWishartInput input = read_input(object);
  const bayests::VarTvpWishartDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarTvpWishartSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VarTvpWishartLogLik)]]
Rcpp::List VarTvpWishartLogLik(Rcpp::List object) {

  const bayests::VarTvpWishartInput input = read_input(object);
  const bayests::VarTvpWishartDraws draws = read_draws_for_loglik(object);

  const arma::mat loglik = bayests::VarTvpWishartSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata, p = 2,
                           deterministic = "const",
                           tvp = TRUE, iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1 / 9, v_i_det = 1 / 100, shape = 3, rate = .0001),
                     sigma = list(df = 3, scale = 1))

object <- add_initial_values(object)

object <- .VarTvpWishartCoefficients(object)
object <- .VarTvpWishartLogLik(object)

*/
