#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_tvp_ald.h"

// VAR estimated at a conditional quantile whose coefficients follow a random
// walk: VarNormalAld with the coefficients turned into a path, exactly as
// VarTvpStochvol is VarNormalStochvol with the same change. The numerics are
// the vendored BayesTS core; what is here is the translation between the R
// model object and the core's structs. See src/core/VENDORED.md.
//
// As for VarNormalAld there is no covariance block and no forecast, so there is
// no .VarTvpAldForecasts beside the two entry points below.

namespace {

using namespace bayests_r;

bayests::VarTvpAldInput read_input(const Rcpp::List &object) {

  bayests::VarTvpAldInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "z", input.train.z);
    }
  }

  const arma::uword tt = input.train.y.n_elem > 0 && input.spec.k > 0
                           ? input.train.periods(input.spec.k)
                           : 0;

  const Rcpp::List initial = has(object, "initial") ? Rcpp::List(object["initial"]) : Rcpp::List();
  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();

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
  if (has(priors, "u_scale")) {
    input.u_scale_prior = read_gamma_prior(Rcpp::List(priors["u_scale"]));
  }

  read_path_if_present(initial, "a", input.train.nparams(), tt, input.initial.a);
  read_mat_if_present(initial, "a_sigma_inv", input.initial.a_sigma_inv);
  read_vec_if_present(initial, "a_init", input.initial.a_init);
  read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);

  read_mat_if_present(initial, "w", input.initial.w);
  read_vec_if_present(initial, "u_scale", input.initial.u_scale);

  return input;
}

/// Every period is scored under its own coefficients, so `a` is the whole path.
/// The asymmetric Laplace density is marginal of the latent scales, so the
/// precision path is read for its column count alone -- iterations() counts it
/// -- and one period is enough for that.
bayests::VarTvpAldDraws read_draws_for_loglik(const Rcpp::List &object,
                                              const bayests::VarTvpAldInput &input) {

  bayests::VarTvpAldDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];

  if (has(posterior, "a")) {
    read_draws_if_present(Rcpp::List(posterior["a"]), "coeffs", draws.a);
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

  return draws;
}

Rcpp::List write_draws(const bayests::VarTvpAldDraws &draws) {

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

  posteriors["u_scale"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_scale));
  posteriors["u_omega_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_omega_inv));
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VarTvpAldCoefficients)]]
Rcpp::List VarTvpAldCoefficients(Rcpp::List object) {

  const bayests::VarTvpAldInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarTvpAldDraws draws =
    bayests::VarTvpAldSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VarTvpAldLogLik)]]
Rcpp::List VarTvpAldLogLik(Rcpp::List object) {

  const bayests::VarTvpAldInput input = read_input(object);
  const bayests::VarTvpAldDraws draws = read_draws_for_loglik(object, input);

  const arma::mat loglik = bayests::VarTvpAldSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvarmodel(data = us_macrodata,
                           p = 1, deterministic = "const",
                           error = "ald", quantile = 0.8, tvp = TRUE,
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1, shape = 3, rate = .0001),
                     sigma = list(shape = 3, rate = .01))

object <- add_initial_values(object)

object <- .VarTvpAldCoefficients(object)
object <- .VarTvpAldLogLik(object)

*/
