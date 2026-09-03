#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_tvp_wishart.h"

// VEC whose coefficients follow a random walk, cointegration vectors included,
// with a Wishart prior on the error precision. It carries no covariance block:
// the error covariance is the Wishart precision alone. The numerics are the
// vendored BayesTS core; what is here is the translation between the R model
// object and the core's structs, plus the entry points
// add_posterior_forecasts() and add_posterior_loglik() dispatch to. See
// src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecTvpWishartInput read_input(const Rcpp::List &object) {

  bayests::VecTvpWishartInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

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
    if (bvs) {
      input.a_varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
    }
  }

  if (has(priors, "beta")) {
    input.beta_prior = read_coint_space_prior_tvp(priors["beta"]);
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

  read_path_if_present(initial, "beta", static_cast<arma::uword>(input.spec.n_beta()), tt,
                       input.initial.beta);
  read_vec_if_present(initial, "beta_init", input.initial.beta_init);

  read_mat_if_present(initial, "u_sigma_inv", input.initial.u_sigma_inv);

  return input;
}

/// Every period is evaluated under its own coefficients, so both paths are read
/// whole. The precision does not move with time in this model, so it is read as
/// it stands.
bayests::VecTvpWishartDraws read_draws_for_loglik(const Rcpp::List &object) {

  bayests::VecTvpWishartDraws draws;

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

/// The forecast starts from where the sample ended, so the coefficients and the
/// cointegration vectors are both sliced to their last in-sample period. Their
/// widths differ and neither is `z.n_cols`: `a` counts off the spec, since a VEC
/// stores the loadings in it, and `beta` off the cointegration dimensions.
bayests::VecTvpWishartDraws read_draws_for_forecast(const Rcpp::List &object,
                                                    const bayests::VecTvpWishartInput &input) {

  bayests::VecTvpWishartDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  const arma::uword tt = input.train.periods(input.spec.k);
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
    read_draws_if_present(Rcpp::List(posterior["u_sigma_inv"]), "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// A member left NULL is how "the model did not have that part" is expressed on
/// the R side, which is the convention the core expresses with an empty matrix.
Rcpp::List write_draws(const bayests::VecTvpWishartDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
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

  // The cointegration path. Without it `a` carries only the loadings, so nothing
  // downstream could reconstruct Pi -- and neither the forecast nor the log
  // likelihood, both of which rebuild the loadings' regressors from it, could be
  // computed at all.
  if (draws.has_beta()) {
    posteriors["beta"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.beta));
  }

  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VecTvpWishartCoefficients)]]
Rcpp::List VecTvpWishartCoefficients(Rcpp::List object) {

  const bayests::VecTvpWishartInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming the
  // first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecTvpWishartDraws draws =
    bayests::VecTvpWishartSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VecTvpWishartForecasts)]]
Rcpp::List VecTvpWishartForecasts(Rcpp::List object) {

  const bayests::VecTvpWishartInput input = read_input(object);
  const bayests::VecTvpWishartDraws draws = read_draws_for_forecast(object, input);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecTvpWishartSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VecTvpWishartLogLik)]]
Rcpp::List VecTvpWishartLogLik(Rcpp::List object) {

  const bayests::VecTvpWishartInput input = read_input(object);
  const bayests::VecTvpWishartDraws draws = read_draws_for_loglik(object);

  // Draws by periods already, which is the orientation R wants and the one WAIC
  // and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecTvpWishartSampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

object <- create_bvecmodel(data = us_macrodata,
                           p = 2, const = "unrestricted",
                           r = 1, tvp = TRUE,
                           error = "wishart",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1, shape = 3, rate = .0001),
                     coint = list(rho = .999),
                     sigma = list(df = 3, scale = 1))

object <- add_initial_values(object)

object <- .VecTvpWishartCoefficients(object)
object <- .VecTvpWishartLogLik(object)

object <- add_forecast_input(object, n_ahead = 10)
object <- .VecTvpWishartForecasts(object)

*/
