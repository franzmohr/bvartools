#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_tvp_discount.h"

// VEC whose loadings follow a random walk and whose error covariance drifts,
// estimated in closed form rather than sampled. Everything the comment at the
// head of src/VarTvpDiscount.cpp says about the shape of the posterior, the RNG
// and the log-likelihood stage holds here; what this file adds is the two
// things a VEC needs.
//
// The design the filter runs against is
//
//     x_t = [ beta' w_t , the compact regressors of data$train$x ],
//
// so the `rank` error correction columns go in front of the short-run blocks,
// which is the layout every VEC here stores `a` in and the reason a column of
// this posterior converts to the level VAR through vec_to_var() without being
// rearranged first.
//
// The cointegration matrix is **given rather than estimated**, and that is what
// makes the model conjugate: beta' w_t is a parameter times data, so a VEC that
// draws beta has a design that moves with the draw and no closed form at all.
// Hold beta and the design is data, every equation shares it, and the whole of
// the discounted VAR applies unchanged. It is read from initial$beta, where
// add_initial_values() puts the space the run conditions on, and written back
// out beside the loadings, which mean nothing without it: alpha and beta are
// identified only as their product.

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

bayests::VecTvpDiscountInput read_input(const Rcpp::List &object) {

  bayests::VecTvpDiscountInput input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "w", input.train.w);
      // The short-run blocks alone -- lagged differences, unmodelled variables
      // and unrestricted deterministic terms -- as for VecKlgs2010. Never `z`.
      read_mat_if_present(train, "x", input.train.x);
    }
    if (has(data, "forecast")) {
      const Rcpp::List forecast = data["forecast"];
      // In the level layout every VEC forecasts in, not the differenced one
      // data$train$x uses.
      read_forecast_regressors(forecast, input.spec.k, input.forecast.x);
    }
    if (has(data, "test")) {
      const Rcpp::List test = data["test"];
      read_mat_if_present(test, "y", input.test.y);
    }
  }

  const Rcpp::List initial = has(object, "initial") ? Rcpp::List(object["initial"]) : Rcpp::List();
  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();

  if (input.use_beta() && has(initial, "beta")) {
    // Stored as vec of a k_beta x rank matrix, the layout every VEC here keeps
    // its space in.
    const arma::vec flat = Rcpp::as<arma::vec>(initial["beta"]);
    const arma::uword k_beta = static_cast<arma::uword>(input.spec.k_beta);
    const arma::uword rank = static_cast<arma::uword>(input.spec.rank);
    if (k_beta > 0 && flat.n_elem == k_beta * rank) {
      input.beta = arma::reshape(flat, k_beta, rank);
    } else {
      // Left in whatever shape the object implies, for validate() to refuse
      // with the message that names both counts. Reshaping a vector of the
      // wrong length here would pad it with zeros and hand the filter a
      // cointegration matrix nobody described.
      input.beta = flat;
    }
  }

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

bayests::VecTvpDiscountPosterior read_posterior(const Rcpp::List &object) {

  bayests::VecTvpDiscountPosterior posterior;

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
  if (has(stored, "beta")) {
    const Rcpp::List beta = stored["beta"];
    if (has(beta, "coeffs")) {
      posterior.beta = arma::vectorise(Rcpp::as<arma::mat>(beta["coeffs"]));
    }
  }

  return posterior;
}

Rcpp::List write_posterior(const bayests::VecTvpDiscountPosterior &posterior) {

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("a") = Rcpp::List::create(
      Rcpp::Named("mean") = draws_to_r(posterior.a),
      Rcpp::Named("scale") = draws_to_r(posterior.a_scale),
      Rcpp::Named("cov") = draws_to_r(posterior.a_cov)),
    Rcpp::Named("u_sigma") = Rcpp::List::create(
      Rcpp::Named("scale") = draws_to_r(posterior.u_sigma)),
    Rcpp::Named("df") = draws_to_r(arma::mat(posterior.df.t())));

  // One row, at the name every other VEC writes its draws of the space to: it
  // did not move and was not estimated, so there is one of it rather than one
  // per period or one per draw.
  if (posterior.has_beta()) {
    out.push_back(Rcpp::List::create(
      Rcpp::Named("coeffs") = draws_to_r(arma::mat(posterior.beta))), "beta");
  }

  return out;
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

arma::uword forecast_draws(const bayests::VarSpec &spec) {
  if (spec.iterations <= 0) {
    Rcpp::stop("A forecast of a discounted model is drawn from its posterior, and 'iterations' "
               "says how many paths to draw, so it must be positive.");
  }
  return static_cast<arma::uword>(spec.iterations);
}

void require_posterior(const bayests::VecTvpDiscountPosterior &posterior) {
  if (posterior.periods() == 0) {
    Rcpp::stop("The model carries no estimated posterior. Call add_posterior_coefficients() "
               "before this step.");
  }
}

} // namespace

// [[Rcpp::export(.VecTvpDiscountCoefficients)]]
Rcpp::List VecTvpDiscountCoefficients(Rcpp::List object) {

  const bayests::VecTvpDiscountInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VecTvpDiscountPosterior posterior =
    bayests::VecTvpDiscountEstimator().estimate(input, reporter);

  return with_posterior(object, write_posterior(posterior));
}

// [[Rcpp::export(.VecTvpDiscountForecasts)]]
Rcpp::List VecTvpDiscountForecasts(Rcpp::List object) {

  const bayests::VecTvpDiscountInput input = read_input(object);
  const bayests::VecTvpDiscountPosterior posterior = read_posterior(object);
  require_posterior(posterior);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecTvpDiscountEstimator().forecast(input, posterior, forecast_draws(input.spec),
                                                reporter);

  return with_forecast_member(object, "forecasts", Rcpp::wrap(draws_to_r(forecast.values)));
}

// [[Rcpp::export(.VecTvpDiscountLogLik)]]
Rcpp::List VecTvpDiscountLogLik(Rcpp::List object) {

  const bayests::VecTvpDiscountInput input = read_input(object);

  bvartools::RcppReporter reporter;

  const bayests::VecTvpDiscountPosterior posterior =
    bayests::VecTvpDiscountEstimator().estimate(input, reporter);

  // The density of data$train$y, which for a VEC is in differences. That is the
  // same number as the density of the level it implies -- given the past,
  // y_t = y_{t-1} + dy_t has unit Jacobian -- so it compares with a VAR in
  // levels on the same sample as well as with another VEC.
  const arma::mat loglik = bayests::VecTvpDiscountEstimator().log_likelihood(input, posterior);

  return with_posterior_element(object, "loglik", Rcpp::wrap(loglik));
}

// [[Rcpp::export(.VecTvpDiscountScore)]]
Rcpp::List VecTvpDiscountScore(Rcpp::List object) {

  const bayests::VecTvpDiscountInput input = read_input(object);
  const bayests::VecTvpDiscountPosterior posterior = read_posterior(object);
  require_posterior(posterior);

  // Drawn, where the discounted VAR's score is exact, and the reason is the
  // layout rather than the model: carrying the filter through the realised
  // values would need the error correction term of each scored period, and a
  // VEC's forecast regressors are levels, from which w_t cannot be recovered.
  // So this is draws x scored periods like the three sampling VECs', and
  // 'iterations' says how many.
  const arma::mat score =
    bayests::VecTvpDiscountEstimator().predictive_log_density(input, posterior,
                                                              forecast_draws(input.spec));

  return with_forecast_member(object, "loglik", Rcpp::wrap(score));
}

/*** R

data("e6")

object <- create_bvecmodel(data = e6 * 100, p = 1, r = 1, const = "unrestricted",
                           algorithm = "discount", delta_beta = 0.98,
                           iterations = 20, burnin = 0, thin = 1)

object <- add_priors(object,
                     coef = list(v_i = 1 / 9, v_i_det = 1 / 100),
                     sigma = list(df = "k", scale = 1))

object <- add_initial_values(object)

object <- .VecTvpDiscountCoefficients(object)
object <- .VecTvpDiscountLogLik(object)

*/
