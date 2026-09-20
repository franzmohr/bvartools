#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/var_tvp_discount.h"

#include <stdexcept>

// VAR whose coefficients follow a random walk and whose error covariance
// drifts, estimated in closed form rather than sampled: the matrix normal
// dynamic linear model of West & Harrison (1997, ch. 16) with the discounted
// Wishart of Uhlig (1997). The numerics are the vendored BayesTS core; what is
// here is the translation between the R model object and the core's structs.
// See src/core/VENDORED.md.
//
// Three things differ from the eight samplers beside it, and all three follow
// from this model having an answer rather than a chain. They are the same three
// that differ in the BayesTS command line, so a model estimated here and one
// estimated on the file come out identical.
//
//   * The posterior is one column per period of a closed form, not one row per
//     draw, so it is written under `mean`, `scale`, `cov`, `scale` and `df`
//     rather than under `coeffs`, and carries no mcpar. add_posterior_coefficients()
//     leaves it alone for exactly that reason.
//   * Nothing in the estimation consumes the RNG, so two runs agree to the bit
//     whatever the seed is. A forecast does consume it, being i.i.d. draws from
//     the closed form, which is what the seed is still there for.
//   * `.VarTvpDiscountLogLik` re-runs the filter rather than reading the stored
//     posterior back. The pointwise score is what the filter produces on its way
//     through and it is not part of what the posterior blocks hold; re-running
//     costs one deterministic pass over the sample and is the same arithmetic on
//     the same input, not a second estimate of it.

namespace {

using namespace bayests_r;

/// The inverse of draws_to_r() for a posterior path: R stores one row per
/// period, the core wants one column per period.
arma::mat path_from_r(const Rcpp::List &list, const char *name) {
  arma::mat out;
  if (has(list, name)) {
    out = arma::trans(Rcpp::as<arma::mat>(list[name]));
  }
  return out;
}

bayests::VarTvpDiscountInput read_input(const Rcpp::List &object) {

  bayests::VarTvpDiscountInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      // The compact regressors, never the SUR matrix `z`: every equation of
      // this model shares its regressors, which is what makes the posterior
      // conjugate, and validate() refuses a file that carries `z` alone.
      read_mat_if_present(train, "x", input.train.x);
    }
    if (has(data, "forecast")) {
      const Rcpp::List forecast = data["forecast"];
      read_forecast_regressors(forecast, input.spec.k, input.forecast.x);
    }
    if (has(data, "test")) {
      const Rcpp::List test = data["test"];
      read_mat_if_present(test, "y", input.test.y);
    }
  }

  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();

  if (input.use_a() && has(priors, "a")) {
    input.a_prior = read_matrix_normal_prior(Rcpp::List(priors["a"]));
  }

  if (has(priors, "u_sigma")) {
    const Rcpp::List prior_u_sigma = priors["u_sigma"];
    input.u_sigma_prior.df = optional_int(prior_u_sigma, "df", 0);
    read_mat_if_present(prior_u_sigma, "scale", input.u_sigma_prior.scale);
  }

  return input;
}

/// The stored posterior, in the shape the core's own file reader rebuilds it.
///
/// `loglik` is not among these and cannot be: it is not one of the blocks the
/// posterior is written as, which is why the log-likelihood entry point re-runs
/// the filter instead of coming through here. Forecasting and scoring need only
/// the last period of `a`, `a_cov`, `u_sigma` and `df`, all four of which are.
bayests::VarTvpDiscountPosterior read_posterior(const Rcpp::List &object) {

  bayests::VarTvpDiscountPosterior posterior;

  if (!has(object, "posterior")) {
    return posterior;
  }

  const Rcpp::List stored = object["posterior"];

  if (has(stored, "a")) {
    const Rcpp::List a = stored["a"];
    posterior.a = path_from_r(a, "mean");
    posterior.a_scale = path_from_r(a, "scale");
    posterior.a_cov = path_from_r(a, "cov");
  }
  if (has(stored, "u_sigma")) {
    posterior.u_sigma = path_from_r(Rcpp::List(stored["u_sigma"]), "scale");
  }
  if (has(stored, "df")) {
    posterior.df = arma::vectorise(Rcpp::as<arma::mat>(stored["df"]));
  }

  return posterior;
}

Rcpp::List write_posterior(const bayests::VarTvpDiscountPosterior &posterior) {

  return Rcpp::List::create(
    Rcpp::Named("a") = Rcpp::List::create(
      Rcpp::Named("mean") = draws_to_r(posterior.a),
      Rcpp::Named("scale") = draws_to_r(posterior.a_scale),
      // The regressor side of the coefficient covariance, n_reg squared per
      // period, kept whole rather than reduced to `scale`: the two together are
      // the joint posterior, and a forecast needs the joint and not the
      // marginal bands.
      Rcpp::Named("cov") = draws_to_r(posterior.a_cov)),
    // A covariance, which is why it is not `u_sigma_inv` where every sampler
    // here puts a precision.
    Rcpp::Named("u_sigma") = Rcpp::List::create(
      Rcpp::Named("scale") = draws_to_r(posterior.u_sigma)),
    Rcpp::Named("df") = draws_to_r(arma::mat(posterior.df.t())));
}

/// The model with its posterior replaced, everything else left as it was.
///
/// Not the four named elements the samplers' entry points rebuild: a discounted
/// model need not have an `initial` at all -- nothing iterates, so there is
/// nowhere for a starting value to be the start of, and only a VEC of positive
/// rank puts the cointegration matrix there -- and naming an element that is
/// absent is an error rather than a NULL.
Rcpp::List with_posterior(const Rcpp::List &object, const Rcpp::List &posterior) {
  Rcpp::List result(Rf_shallow_duplicate(object));
  if (result.containsElementNamed("posterior")) {
    result["posterior"] = posterior;
  } else {
    result.push_back(posterior, "posterior");
  }
  return result;
}

/// How many i.i.d. draws a stage that needs them takes. `iterations` keeps its
/// name and loses its chain: there is no sweep, so the number is simply how many
/// draws of the posterior to take.
arma::uword forecast_draws(const bayests::VarSpec &spec) {
  if (spec.iterations <= 0) {
    Rcpp::stop("A forecast of a discounted model is drawn from its posterior, and 'iterations' "
               "says how many paths to draw, so it must be positive.");
  }
  return static_cast<arma::uword>(spec.iterations);
}

void require_posterior(const bayests::VarTvpDiscountPosterior &posterior) {
  if (posterior.periods() == 0) {
    Rcpp::stop("The model carries no estimated posterior. Call add_posterior_coefficients() "
               "before this step.");
  }
}

} // namespace

// [[Rcpp::export(.VarTvpDiscountCoefficients)]]
Rcpp::List VarTvpDiscountCoefficients(Rcpp::List object) {

  const bayests::VarTvpDiscountInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarTvpDiscountPosterior posterior =
    bayests::VarTvpDiscountEstimator().estimate(input, reporter);

  return with_posterior(object, write_posterior(posterior));
}

// [[Rcpp::export(.VarTvpDiscountForecasts)]]
Rcpp::List VarTvpDiscountForecasts(Rcpp::List object) {

  const bayests::VarTvpDiscountInput input = read_input(object);
  const bayests::VarTvpDiscountPosterior posterior = read_posterior(object);
  require_posterior(posterior);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VarTvpDiscountEstimator().forecast(input, posterior, forecast_draws(input.spec),
                                                reporter);

  return with_forecast_member(object, "forecasts", Rcpp::wrap(draws_to_r(forecast.values)));
}

// [[Rcpp::export(.VarTvpDiscountLogLik)]]
Rcpp::List VarTvpDiscountLogLik(Rcpp::List object) {

  const bayests::VarTvpDiscountInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VarTvpDiscountPosterior posterior =
    bayests::VarTvpDiscountEstimator().estimate(input, reporter);

  // One row rather than one per draw: the parameters are integrated out
  // exactly, so the sum of this is the log marginal likelihood of the sample
  // and not an estimate of it. selection_criteria() reports it as LML.
  const arma::mat loglik = bayests::VarTvpDiscountEstimator().log_likelihood(input, posterior);

  return with_posterior_element(object, "loglik", Rcpp::wrap(loglik));
}

// [[Rcpp::export(.VarTvpDiscountScore)]]
Rcpp::List VarTvpDiscountScore(Rcpp::List object) {

  const bayests::VarTvpDiscountInput input = read_input(object);
  const bayests::VarTvpDiscountPosterior posterior = read_posterior(object);
  require_posterior(posterior);

  // One row again, and for the same reason: the filter carries itself through
  // the realised values, each scored period conditioning on the ones before it,
  // so there is nothing to average over.
  const arma::mat score =
    bayests::VarTvpDiscountEstimator().predictive_log_density(input, posterior);

  return with_forecast_member(object, "loglik", Rcpp::wrap(score));
}

/*** R

data("at_macrodata")

object <- create_bvarmodel(data = at_macrodata[["domestic"]][, c("y", "Dp")] * 100,
                           p = 1, deterministic = "const", algorithm = "discount",
                           delta_beta = 0.99, delta_sigma = 0.98,
                           iterations = 20, burnin = 0, thin = 1)

object <- add_priors(object,
                     coef = list(v_i = 1 / 9, v_i_det = 1 / 100),
                     sigma = list(df = "k", scale = 1))

object <- add_initial_values(object)

object <- .VarTvpDiscountCoefficients(object)
object <- .VarTvpDiscountLogLik(object)

*/
