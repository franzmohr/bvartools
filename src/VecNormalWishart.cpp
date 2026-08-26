#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_normal_wishart.h"

// VEC with a normal prior on non-cointegration coefficients and a prior on the
// cointegration space and a Wishart prior on the error precision.
// The numerics are the vendored BayesTS core; what is here is the
// translation between the R model object and the core's structs, plus the
// entry points add_posterior_forecasts() and add_posterior_loglik()
// dispatch to. See src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecNormalWishartInput read_input(const Rcpp::List &object) {

  bayests::VecNormalWishartInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "w", input.train.w);
      read_mat_if_present(train, "z", input.train.z);
    }
    if (has(data, "forecast")) {
      const Rcpp::List forecast = data["forecast"];
      read_mat_if_present(forecast, "z", input.forecast.z);
    }
  }

  // Only draw_coefficients() needs the priors and the initial values, so a
  // missing one is left for validate() to complain about if it turns out to
  // matter: an object holding nothing but a fitted posterior can still be
  // forecast from.
  if (has(object, "priors")) {
    const Rcpp::List priors = object["priors"];
    
    if (has(priors, "a")) {
      const Rcpp::List prior_a = priors["a"];
      input.a_prior = read_normal_prior(prior_a);
      if (input.spec.uses_varsel()) {
        // R keeps the selection prior alongside the normal one rather than in
        // a group of its own.
        input.varsel_prior = read_varsel_prior(prior_a, input.spec.varsel);
      }
    }
    
    if (has(priors, "beta")) {
      const Rcpp::List prior_beta = priors["beta"];
      input.beta_prior = read_coint_space_prior_constant(prior_beta);
    }
    
    if (has(priors, "u_sigma")) {
      const Rcpp::List prior_u_sigma = priors["u_sigma"];
      input.u_sigma_prior.df = optional_int(prior_u_sigma, "df", 0);
      read_mat_if_present(prior_u_sigma, "scale", input.u_sigma_prior.scale);
    }
  }

  if (has(object, "initial")) {
    const Rcpp::List initial = object["initial"];
    read_vec_if_present(initial, "a", input.initial.a);
    read_vec_if_present(initial, "a_lambda", input.initial.a_lambda);
    read_vec_if_present(initial, "beta", input.initial.beta);
    read_mat_if_present(initial, "u_sigma_inv", input.initial.u_sigma_inv);
  }

  return input;
}

/// Both the forecast and the log likelihood read every draw as it stands:
/// neither the coefficients nor the precision move with time in this model.
bayests::VecNormalWishartDraws read_draws(const Rcpp::List &object) {

  bayests::VecNormalWishartDraws draws;

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
    const Rcpp::List posterior_beta = posterior["beta"];
    read_draws_if_present(posterior_beta, "coeffs", draws.beta);
  }
  if (has(posterior, "u_sigma_inv")) {
    const Rcpp::List posterior_u_sigma_inv = posterior["u_sigma_inv"];
    read_draws_if_present(posterior_u_sigma_inv, "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

/// A member left NULL is how "the model did not have that part" is expressed on
/// the R side, which is the convention the core expresses with an empty matrix.
Rcpp::List write_draws(const bayests::VecNormalWishartDraws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
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
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VecNormalWishartCoefficients)]]
Rcpp::List VecNormalWishartCoefficients(Rcpp::List object) {

  const bayests::VecNormalWishartInput input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming
  // the first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecNormalWishartDraws draws =
    bayests::VecNormalWishartSampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VecNormalWishartForecasts)]]
Rcpp::List VecNormalWishartForecasts(Rcpp::List object) {

  const bayests::VecNormalWishartInput input = read_input(object);
  const bayests::VecNormalWishartDraws draws = read_draws(object);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecNormalWishartSampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VecNormalWishartLogLik)]]
Rcpp::List VecNormalWishartLogLik(Rcpp::List object) {

  const bayests::VecNormalWishartInput input = read_input(object);
  const bayests::VecNormalWishartDraws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one
  // WAIC and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecNormalWishartSampler().log_likelihood(input, draws);

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
                           error = "wishart",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 3, scale = 1))

object <- add_initial_values(object)

object <- .VecNormalWishartCoefficients(object)
object <- .VecNormalWishartLogLik(object)

object <- add_forecast_input_data(object, n_ahead = 10)
object <- .VecNormalWishartForecasts(object)

*/
