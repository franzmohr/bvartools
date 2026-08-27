#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests_reporter.h"
#include "bayests/vec_klgs_2010.h"

// The non-SUR reading of Koop, Leon-Gonzalez and Strachan (2010): the same
// posterior .VecNormalWishartCoefficients() draws, obtained without forming the
// SUR design matrix. It reads data$train$x, the compact regressors
// create_bvecmodel() already assembles, in place of data$train$z; everything
// else about the model object is the same, and no selection scheme is
// available -- the core's validate() rejects one and says why.
//
// The numerics are the vendored BayesTS core; what is here is the translation
// between the R model object and the core's structs, plus the entry points
// add_posterior_forecasts() and add_posterior_loglik() dispatch to. See
// src/core/VENDORED.md.

namespace {

using namespace bayests_r;

bayests::VecKlgs2010Input read_input(const Rcpp::List &object) {

  bayests::VecKlgs2010Input input;

  const Rcpp::List model = object["model"];
  input.spec = read_spec(model);

  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      const Rcpp::List train = data["train"];
      read_mat_if_present(train, "y", input.train.y);
      read_mat_if_present(train, "w", input.train.w);
      // The compact regressors, not the SUR ones. A model object carries both,
      // so reading the wrong name here would not fail -- it would run the
      // sampler on a matrix with tt*k rows and produce numbers.
      read_mat_if_present(train, "x", input.train.x);
    }
    if (has(data, "forecast")) {
      // Still the SUR layout, and still in levels: the forecast is the level
      // VAR's for every VEC, and it is VarNormalWishartSampler that reads this.
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
    read_vec_if_present(initial, "beta", input.initial.beta);
    read_mat_if_present(initial, "u_sigma_inv", input.initial.u_sigma_inv);
  }

  return input;
}

/// Both the forecast and the log likelihood read every draw as it stands:
/// neither the coefficients nor the precision move with time in this model.
bayests::VecKlgs2010Draws read_draws(const Rcpp::List &object) {

  bayests::VecKlgs2010Draws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  if (has(posterior, "a")) {
    const Rcpp::List posterior_a = posterior["a"];
    read_draws_if_present(posterior_a, "coeffs", draws.a);
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
Rcpp::List write_draws(const bayests::VecKlgs2010Draws &draws) {

  Rcpp::List posteriors = Rcpp::List::create(Rcpp::Named("a") = R_NilValue,
                                             Rcpp::Named("beta") = R_NilValue,
                                             Rcpp::Named("u_sigma_inv") = R_NilValue);

  if (draws.has_a()) {
    posteriors["a"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.a));
  }
  if (draws.has_beta()) {
    posteriors["beta"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.beta));
  }
  posteriors["u_sigma_inv"] = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.u_sigma_inv));

  return posteriors;
}

} // namespace

// [[Rcpp::export(.VecKlgs2010Coefficients)]]
Rcpp::List VecKlgs2010Coefficients(Rcpp::List object) {

  const bayests::VecKlgs2010Input input = read_input(object);

  // Throttled Rcpp::checkUserInterrupt(); silent unless asked to report.
  bvartools::RcppReporter reporter;

  // The sampler validates the input and throws std::invalid_argument naming
  // the first inconsistency it finds; Rcpp turns that into an R error.
  const bayests::VecKlgs2010Draws draws =
    bayests::VecKlgs2010Sampler().draw_coefficients(input, reporter);

  return Rcpp::List::create(Rcpp::Named("data") = object["data"],
                            Rcpp::Named("model") = object["model"],
                            Rcpp::Named("initial") = object["initial"],
                            Rcpp::Named("priors") = object["priors"],
                            Rcpp::Named("posterior") = write_draws(draws));
}

// [[Rcpp::export(.VecKlgs2010Forecasts)]]
Rcpp::List VecKlgs2010Forecasts(Rcpp::List object) {

  const bayests::VecKlgs2010Input input = read_input(object);
  const bayests::VecKlgs2010Draws draws = read_draws(object);

  bvartools::RcppReporter reporter;

  const bayests::ForecastDraws forecast =
    bayests::VecKlgs2010Sampler().forecast(input, draws, reporter);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(draws_to_r(forecast.values), "forecast");
  object["posterior"] = posterior;

  return object;
}

// [[Rcpp::export(.VecKlgs2010LogLik)]]
Rcpp::List VecKlgs2010LogLik(Rcpp::List object) {

  const bayests::VecKlgs2010Input input = read_input(object);
  const bayests::VecKlgs2010Draws draws = read_draws(object);

  // Draws by periods already, which is the orientation R wants and the one
  // WAIC and PSIS-LOO expect; no transpose at this boundary.
  const arma::mat loglik = bayests::VecKlgs2010Sampler().log_likelihood(input, draws);

  Rcpp::List posterior = object["posterior"];
  posterior.push_back(loglik, "loglik");
  object["posterior"] = posterior;

  return object;
}

/*** R

data("us_macrodata")

model <- function() {
  object <- create_bvecmodel(data = us_macrodata,
                             p = 2, const = "unrestricted",
                             r = 1,
                             error = "wishart",
                             iterations = 20, burnin = 10)

  object <- add_priors(object,
                       coef = list(v_i = 1),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 1))

  add_initial_values(object)
}

object <- .VecKlgs2010Coefficients(model())
object <- .VecKlgs2010LogLik(object)

## The same posterior as the SUR sampler, from the same seed. This is the check
## worth running after touching either of them.
set.seed(4711); klgs <- .VecKlgs2010Coefficients(model())
set.seed(4711); sur <- .VecNormalWishartCoefficients(model())
max(abs(klgs$posterior$a$coeffs - sur$posterior$a$coeffs))

## .VecKlgs2010Forecasts() is not reachable from a 'bvecmodel' on its own: the
## R side forecasts a VEC by converting it with vec_to_var() first, so the
## forecast regressors this entry point wants -- data$forecast$z, in levels --
## are assembled there rather than here.

*/
