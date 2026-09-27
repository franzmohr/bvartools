#ifndef BVARTOOLS_BAYESTS_R_IO_H
#define BVARTOOLS_BAYESTS_R_IO_H

#include <RcppArmadillo.h>

#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/results.h"
#include "bayests/spec.h"

#include <cmath>
#include <string>
#include <type_traits>
#include <utility>

// Translation between the R model object and the structs the vendored BayesTS
// core takes. This is the R counterpart of that project's src/io/hdf5/ layer:
// the only place that knows both an `Rcpp::List` and a `bayests::` struct, so
// that the samplers know neither.
//
// The pieces every model shares live here; a model's own reader lives next to
// its binding. Two conventions are worth stating once, because getting either
// wrong is silent rather than loud:
//
//   * Draws run along the columns inside the core and along the rows in R.
//     Everything crossing this boundary goes through draws_to_r() or
//     draws_from_r().
//   * Coefficient positions are counted from one in R and from zero in the
//     core. read_varsel_prior() is the only place that converts.
namespace bayests_r
{

inline bool has(const Rcpp::List &list, const char *name)
{
  return list.containsElementNamed(name) && !Rf_isNull(list[name]);
}

inline int optional_int(const Rcpp::List &list, const char *name, int fallback)
{
  return has(list, name) ? Rcpp::as<int>(list[name]) : fallback;
}

inline bool optional_bool(const Rcpp::List &list, const char *name, bool fallback)
{
  return has(list, name) ? Rcpp::as<bool>(list[name]) : fallback;
}

inline std::string optional_string(const Rcpp::List &list, const char *name,
                                   const std::string &fallback)
{
  return has(list, name) ? Rcpp::as<std::string>(list[name]) : fallback;
}

inline void read_double_if_present(const Rcpp::List &list, const char *name, double &out)
{
  if (has(list, name)) {
    out = Rcpp::as<double>(list[name]);
  }
}

inline void read_mat_if_present(const Rcpp::List &list, const char *name, arma::mat &out)
{
  if (has(list, name)) {
    out = Rcpp::as<arma::mat>(list[name]);
  }
}

inline void read_vec_if_present(const Rcpp::List &list, const char *name, arma::vec &out)
{
  if (has(list, name)) {
    out = Rcpp::as<arma::vec>(list[name]);
  }
}

/// The out-of-sample regressors, in the compact layout ForecastData::x is
/// written in: one row per horizon, one column per regressor.
///
/// Takes `x` when the list has it. An object built by an older version of this
/// package carries `z` instead -- the same regressors kroneckered up with I_k,
/// at k times the rows and k times the columns -- and is compacted back on the
/// way in, so a fitted model saved to disk before the layout changed still
/// forecasts. `k` is what decides which of the two a `z` is, so it comes from
/// the model specification rather than from the matrix's own shape.
///
/// Leaves `out` alone when the list has neither, which is what a model with no
/// forecast requested looks like; require_forecast_regressors() in the core is
/// what turns that into an error for a model that needed them.
inline void read_forecast_regressors(const Rcpp::List &list, const int k, arma::mat &out)
{
  if (has(list, "x")) {
    out = Rcpp::as<arma::mat>(list["x"]);
    return;
  }
  if (!has(list, "z")) {
    return;
  }

  const arma::mat sur = Rcpp::as<arma::mat>(list["z"]);
  if (k <= 0) {
    Rcpp::stop("data$forecast$z is the SUR layout and can only be read against a positive k");
  }
  const arma::uword width = static_cast<arma::uword>(k);
  if (sur.n_rows % width != 0 || sur.n_cols % width != 0) {
    Rcpp::stop("data$forecast$z is the SUR layout, so both of its dimensions have to be "
               "multiples of k");
  }

  // z is kron(x, I_k), so the block at row i and column j is x(i, j) * I_k and
  // x(i, j) is the one element of it that every block has in the same place.
  out.set_size(sur.n_rows / width, sur.n_cols / width);
  for (arma::uword i = 0; i < out.n_rows; i++) {
    for (arma::uword j = 0; j < out.n_cols; j++) {
      out(i, j) = sur(i * width, j * width);
    }
  }
}

/// Draws as the core keeps them, from the row-per-draw matrix R keeps.
inline void read_draws_if_present(const Rcpp::List &list, const char *name, arma::mat &out)
{
  if (has(list, name)) {
    out = arma::trans(Rcpp::as<arma::mat>(list[name]));
  }
}

/// Draws as R expects them, from the column-per-draw matrix the core returns.
inline arma::mat draws_to_r(const arma::mat &draws)
{
  return arma::trans(draws);
}

/// A block's draws with what the non-centred parameterisation adds beside its
/// `sigma`: the signed standard deviation and the log ordinates at zero the
/// Savage-Dickey test for time variation is built from. Unchanged for a
/// centred block. Shared by the bindings of the models that offer `omega_v`.
inline Rcpp::List with_noncentred(Rcpp::List block, const bayests::NoncentredStateDraws &nc)
{
  if (!nc.empty()) {
    block["omega"] = draws_to_r(nc.omega);
    block["omega_log_zero"] = draws_to_r(nc.log_zero);
    block["omega_log_zero_joint"] = draws_to_r(nc.log_zero_joint);
  }
  return block;
}

/// A time-varying starting value. R stores the whole path flat; the core wants
/// it `width` by `tt`, one column per period.
inline void read_path_if_present(const Rcpp::List &list, const char *name,
                                 arma::uword width, arma::uword tt, arma::mat &out)
{
  if (has(list, name)) {
    out = arma::reshape(Rcpp::as<arma::vec>(list[name]), width, tt);
  }
}

/// The last period of a stored path, for the samplers that forecast from where
/// the sample ended: columns [(tt - 1) * width, tt * width) of the row-per-draw
/// matrix, handed over as one column per draw.
inline void read_draws_last_period_if_present(const Rcpp::List &list, const char *name,
                                              arma::uword tt, arma::uword width, arma::mat &out)
{
  if (!has(list, name)) {
    return;
  }
  const arma::mat stored = Rcpp::as<arma::mat>(list[name]);
  if (tt == 0 || width == 0 || stored.n_cols < tt * width) {
    Rcpp::stop("'%s' holds %d columns, too few for the last of %d periods",
               name, static_cast<int>(stored.n_cols), static_cast<int>(tt));
  }
  out = arma::trans(stored.cols((tt - 1) * width, tt * width - 1));
}

/// The model specification. `covar_error` is the value of `model$error` that
/// means "with a covariance block" for the calling model, or `nullptr` for a
/// model that has none.
inline bayests::VarSpec read_spec(const Rcpp::List &model, const char *covar_error = nullptr)
{
  bayests::VarSpec spec;

  spec.k = Rcpp::as<int>(model["k"]);
  spec.iterations = Rcpp::as<int>(model["iterations"]);
  spec.burnin = Rcpp::as<int>(model["burnin"]);
  // Absent unless create_*model() was asked to thin, and 1 then, which keeps
  // every draw after the burn-in.
  spec.thin = optional_int(model, "thin", 1);
  spec.p = optional_int(model, "p", 0);
  spec.m = optional_int(model, "m", 0);
  spec.s = optional_int(model, "s", 0);
  spec.n = optional_int(model, "n", 0);
  // Endogenous variables, ordered first, whose equations carry no coefficients
  // at all. Absent from every model but one create_bvarmodel() was given `iid`,
  // and refused by the sampler's own validate() where no sampler reads it, so
  // it is read for every algorithm rather than only the four that do.
  spec.n_iid = optional_int(model, "n_iid", 0);
  spec.n_restricted = optional_int(model, "n_restricted", 0);
  spec.rank = optional_int(model, "rank", 0);
  spec.k_beta = optional_int(model, "k_beta", 0);
  // Unobserved factors. Zero for every model this package has -- a dynamic
  // factor model is dfmtools' -- and read anyway, because the vendored core is
  // a whole mirror of upstream and its VarSpec carries the field.
  spec.n_factors = optional_int(model, "n_factors", 0);
  // Absent until add_forecast_input() has been called.
  spec.h = optional_int(model, "h", 0);
  // The quantile an asymmetric Laplace model estimates. Every other model
  // ignores it, and VarSpec defaults it to 0.5 -- which is why it has to be
  // read rather than left alone: a model asked for the 0.8 quantile whose
  // spec never reaches the sampler estimates the median and says nothing.
  read_double_if_present(model, "quantile", spec.quantile);
  // The two discount factors of VarTvpDiscount and VecTvpDiscount. One apiece
  // is VarSpec's default and the model in which the quantity they govern does
  // not move, so a specification that names neither is a constant coefficient
  // model rather than an unspecified one. Every other algorithm ignores them.
  read_double_if_present(model, "delta_beta", spec.delta_beta);
  read_double_if_present(model, "delta_sigma", spec.delta_sigma);
  spec.varsel = bayests::var_selection_from_string(optional_string(model, "varsel", "none"));
  spec.structural = optional_bool(model, "structural", false);
  // What the core added after this reader was written: the adaptive prior on
  // `a`, the stationarity condition and the steady-state prior. Absent unless
  // add_prior_options() set them, and refused by the sampler's own validate()
  // where no sampler reads them.
  spec.shrinkage = bayests::shrinkage_from_string(optional_string(model, "shrinkage", "none"));
  spec.stationary = optional_bool(model, "stationary", false);
  spec.steady_state = optional_bool(model, "steady_state", false);
  // Whether a time-varying model's forecast carries its random walks over the
  // horizon or holds them at the last sample period. Absent unless
  // add_posterior_forecasts() was given one, and `simulate` then, which is the
  // core's own default.
  spec.forecast_states = bayests::forecast_states_from_string(
    optional_string(model, "forecast_states", "simulate"));

  if (covar_error != nullptr) {
    spec.covar = optional_string(model, "error", "") == covar_error;
  }

  return spec;
}

inline bayests::NormalPrior read_normal_prior(const Rcpp::List &group)
{
  bayests::NormalPrior prior;
  read_vec_if_present(group, "mu", prior.mu);
  read_mat_if_present(group, "v_inv", prior.v_inv);
  return prior;
}

/// The matrix normal prior of a discounted model: the n_design x k coefficient
/// matrix and the n_design square regressor side of its covariance.
///
/// Deliberately not read from `mu` and `v_inv`, which is where a sampler keeps
/// its normal prior over the vectorised coefficients. The two are different
/// objects -- `cov` is a covariance in units of the error covariance, since the
/// full prior covariance is Sigma kronecker cov -- and a model that silently
/// took one for the other would run against a prior nobody stated.
inline bayests::MatrixNormalPrior read_matrix_normal_prior(const Rcpp::List &group)
{
  bayests::MatrixNormalPrior prior;
  read_mat_if_present(group, "mean", prior.mean);
  read_mat_if_present(group, "cov", prior.cov);
  return prior;
}

inline bayests::ConstantCointSpacePrior read_coint_space_prior_constant(const Rcpp::List &group)
{
  bayests::ConstantCointSpacePrior prior;
  read_double_if_present(group, "v_inv", prior.v_inv);
  read_mat_if_present(group, "p_tau_inv", prior.p_tau_inv);
  // G^-1, which only VecNormalStochvol reads; validate() refuses one anywhere
  // else, as cointspace_prior() already does. Absent, the core fixes G from
  // the starting volatilities.
  read_mat_if_present(group, "g_inv", prior.g_inv);
  return prior;
}

/// The cointegration space prior of a model whose cointegration vectors move.
/// `rho` is the autoregression of their state equation.
///
/// Without `rho_min` and `rho_max` it is a fixed hyperparameter, and a group
/// that omits it keeps the struct's default. With them it is drawn, under a
/// uniform prior on that interval, and `rho` becomes the value the chain starts
/// at. Both ends or neither: one alone would leave the sampler to invent the
/// other, and which end is missing changes the model rather than a detail of it.
inline bayests::TvpCointSpacePrior read_coint_space_prior_tvp(const Rcpp::List &group)
{
  bayests::TvpCointSpacePrior prior;
  prior.initial_state = read_normal_prior(group);
  read_double_if_present(group, "rho", prior.rho);

  const bool has_min = has(group, "rho_min");
  const bool has_max = has(group, "rho_max");

  if (has_min != has_max) {
    Rcpp::stop("the prior support of rho needs both ends: priors$beta$%s is missing. "
               "Leave both out to hold rho fixed at priors$beta$rho",
               has_min ? "rho_max" : "rho_min");
  }

  if (has_min) {
    prior.rho_prior.draw = true;
    prior.rho_prior.min = Rcpp::as<double>(group["rho_min"]);
    prior.rho_prior.max = Rcpp::as<double>(group["rho_max"]);
  }

  // The transition of the state equation with rho taken out -- Koop et al.'s
  // informative marginal prior. Absent is the identity.
  read_mat_if_present(group, "p_tau", prior.p_tau);

  return prior;
}

inline bayests::GammaPrior read_gamma_prior(const Rcpp::List &group)
{
  bayests::GammaPrior prior;
  read_vec_if_present(group, "shape", prior.shape);
  read_vec_if_present(group, "rate", prior.rate);
  return prior;
}

/// Positions of the coefficients selection applies to, converted from R's
/// one-based counting to the core's zero-based. Checked before the
/// subtraction: a zero would wrap to an index no bounds check recognises.
inline arma::uvec read_positions(const Rcpp::List &group, const char *name)
{
  const arma::vec one_based = Rcpp::as<arma::vec>(group[name]);

  if (one_based.n_elem > 0 && one_based.min() < 1.0) {
    Rcpp::stop("'%s' holds a position below 1; coefficient positions are counted from one", name);
  }

  return arma::conv_to<arma::uvec>::from(one_based - 1);
}

inline bayests::VarSelPrior read_varsel_prior(const Rcpp::List &group, bayests::VarSelection scheme)
{
  bayests::VarSelPrior prior;
  read_vec_if_present(group, "inprior", prior.inprior);
  if (has(group, "include")) {
    prior.include = read_positions(group, "include");
  }

  if (scheme == bayests::VarSelection::ssvs) {
    read_vec_if_present(group, "tau0", prior.ssvs.tau0);
    read_vec_if_present(group, "tau1", prior.ssvs.tau1);
  }

  return prior;
}

// A copy of `object` whose posterior holds `value` under `name`, replacing an
// element of that name rather than adding a second one beside it.
//
// An Rcpp::List argument wraps the caller's R object, not a copy of it, so
// writing into it changes the object the caller still holds. Only the two lists
// that change are duplicated, and shallowly: the draws already in the posterior
// are shared with the caller's object, not copied.
inline Rcpp::List with_posterior_element(const Rcpp::List &object, const char *name,
                                         const Rcpp::RObject &value)
{
  Rcpp::List result(Rf_shallow_duplicate(object));
  const Rcpp::List current = object["posterior"];
  Rcpp::List posterior(Rf_shallow_duplicate(current));
  if (posterior.containsElementNamed(name)) {
    posterior[name] = value;
  } else {
    posterior.push_back(value, name);
  }
  result["posterior"] = posterior;
  return result;
}

/// One member of the forecast group of the posterior, at
/// posterior$forecast$<name>.
///
/// One level in rather than at posterior$forecast itself, because the group is
/// the place for everything the forecast periods produce: the paths in
/// `forecasts`, the `errors` add_forecast_errors() takes against a test sample,
/// and the `loglik` add_predictive_loglik() scores them with. That is the layout
/// of the model file, where `bayests` writes /posterior/forecast/forecasts, and
/// the leaves are named after what they hold rather than one of them being
/// called `draws` because every member of the group is draws.
///
/// An object fitted before the group existed carries a matrix at
/// posterior$forecast rather than a list. The matrix is the paths, so it is
/// dropped rather than merged into: whatever is about to be written belongs to
/// the draws the caller has just produced, and an `errors` or a `loglik` taken
/// against the old paths would not belong with them anyway.
inline Rcpp::List with_forecast_member(const Rcpp::List &object, const char *name,
                                       const Rcpp::RObject &value)
{
  Rcpp::List result(Rf_shallow_duplicate(object));
  const Rcpp::List current = object["posterior"];
  Rcpp::List posterior(Rf_shallow_duplicate(current));

  Rcpp::List forecast;
  const bool had_group = posterior.containsElementNamed("forecast") &&
                         Rf_isNewList(posterior["forecast"]);
  if (had_group) {
    forecast = Rcpp::List(Rf_shallow_duplicate(posterior["forecast"]));
  }

  if (forecast.containsElementNamed(name)) {
    forecast[name] = value;
  } else {
    forecast.push_back(value, name);
  }

  if (posterior.containsElementNamed("forecast")) {
    posterior["forecast"] = forecast;
  } else {
    posterior.push_back(forecast, "forecast");
  }
  result["posterior"] = posterior;
  return result;
}

// ---------------------------------------------------------------------------
// A panel not observed whole, and the prior options of the constant VARs
// ---------------------------------------------------------------------------
//
// What the core added after the bindings were written, read and written for
// every model in one place: the three constraint groups (data$train, data$test
// and data$forecast, each a list of value, group, row, period, variable and
// weight), the soft groups' prior and start, the adaptive prior on `a`, and
// the steady-state prior on the mean. Each binding calls read_extensions(),
// read_draw_extensions() and write_draw_extensions() whatever its model, and
// what the model's structs have no member for is skipped at compile time. A
// model that does not read a group it is handed refuses it in its validate(),
// as a model file would be refused -- the point of reading it everywhere.

template <typename T, typename = void>
struct has_forecast_constraints : std::false_type {};
template <typename T>
struct has_forecast_constraints<T, std::void_t<decltype(std::declval<T &>().forecast.constraints)>>
    : std::true_type {};

template <typename T, typename = void>
struct has_constraints_prior : std::false_type {};
template <typename T>
struct has_constraints_prior<T, std::void_t<decltype(std::declval<T &>().constraints_prior)>>
    : std::true_type {};

template <typename T, typename = void>
struct has_shrinkage_prior : std::false_type {};
template <typename T>
struct has_shrinkage_prior<T, std::void_t<decltype(std::declval<T &>().a_shrinkage_prior)>>
    : std::true_type {};

template <typename T, typename = void>
struct has_mu_prior : std::false_type {};
template <typename T>
struct has_mu_prior<T, std::void_t<decltype(std::declval<T &>().mu_prior)>> : std::true_type {};

template <typename T, typename = void>
struct has_completion_draws : std::false_type {};
template <typename T>
struct has_completion_draws<T, std::void_t<decltype(std::declval<T &>().constraints_inv)>>
    : std::true_type {};

template <typename T, typename = void>
struct has_shrinkage_draws : std::false_type {};
template <typename T>
struct has_shrinkage_draws<T, std::void_t<decltype(std::declval<T &>().a_shrinkage)>>
    : std::true_type {};

template <typename T, typename = void>
struct has_mu_draws : std::false_type {};
template <typename T>
struct has_mu_draws<T, std::void_t<decltype(std::declval<T &>().mu)>> : std::true_type {};

/// Whole numbers no smaller than `minimum`, shifted down by it, from a numeric
/// R vector: what the file's positions are, counted from one, and its groups,
/// from zero.
inline arma::uvec read_counted(const Rcpp::List &list, const char *name, const double minimum,
                               const char *where)
{
  const arma::vec values = Rcpp::as<arma::vec>(list[name]);
  for (const double v : values) {
    if (!std::isfinite(v) || v != std::floor(v) || v < minimum) {
      Rcpp::stop("%s$%s holds %f, which is not a whole number of at least %d", where, name, v,
                 static_cast<int>(minimum));
    }
  }
  return arma::conv_to<arma::uvec>::from(values - minimum);
}

/// A constraint set from its R list, the six elements of the model file's
/// group; empty when the element is absent or NULL.
inline bayests::Constraints read_constraints(const Rcpp::List &parent, const char *name,
                                             const char *where)
{
  bayests::Constraints c;
  if (!has(parent, name) || Rf_isNull(parent[name])) {
    return c;
  }
  const Rcpp::List group = parent[name];
  for (const char *member : {"value", "group", "row", "period", "variable", "weight"}) {
    if (!has(group, member)) {
      Rcpp::stop("%s has no '%s': a constraint set is value, group, row, period, variable and "
                 "weight together", where, member);
    }
  }
  c.value = Rcpp::as<arma::vec>(group["value"]);
  c.group = read_counted(group, "group", 0.0, where);
  c.row = read_counted(group, "row", 1.0, where);
  c.period = read_counted(group, "period", 1.0, where);
  c.variable = read_counted(group, "variable", 1.0, where);
  c.weight = Rcpp::as<arma::vec>(group["weight"]);
  return c;
}

template <typename Input>
void read_extensions(const Rcpp::List &object, Input &input)
{
  if (has(object, "data")) {
    const Rcpp::List data = object["data"];
    if (has(data, "train")) {
      input.train.constraints = read_constraints(data["train"], "constraints",
                                                 "data$train$constraints");
    }
    if (has(data, "test")) {
      input.test.constraints = read_constraints(data["test"], "constraints",
                                                "data$test$constraints");
    }
    if constexpr (has_forecast_constraints<Input>::value) {
      if (has(data, "forecast")) {
        input.forecast.constraints = read_constraints(data["forecast"], "constraints",
                                                      "data$forecast$constraints");
      }
    }
  }

  const Rcpp::List priors = has(object, "priors") ? Rcpp::List(object["priors"]) : Rcpp::List();
  const Rcpp::List initial = has(object, "initial") ? Rcpp::List(object["initial"]) : Rcpp::List();

  if constexpr (has_constraints_prior<Input>::value) {
    if (has(priors, "constraints")) {
      input.constraints_prior = read_gamma_prior(priors["constraints"]);
    }
    read_vec_if_present(initial, "constraints_inv", input.initial.constraints_inv);
  }
  if constexpr (has_shrinkage_prior<Input>::value) {
    if (has(priors, "a")) {
      const Rcpp::List prior_a = priors["a"];
      if (has(prior_a, "shrinkage")) {
        const Rcpp::List shrinkage = prior_a["shrinkage"];
        if (has(shrinkage, "group")) {
          input.a_shrinkage_prior.group = read_counted(shrinkage, "group", 0.0,
                                                       "priors$a$shrinkage");
        }
        read_vec_if_present(shrinkage, "shape", input.a_shrinkage_prior.shape);
        read_vec_if_present(shrinkage, "rate", input.a_shrinkage_prior.rate);
      }
    }
    read_vec_if_present(initial, "a_shrinkage", input.initial.a_shrinkage);
    read_vec_if_present(initial, "a_local", input.initial.a_local);
  }
  if constexpr (has_mu_prior<Input>::value) {
    if (has(priors, "mu")) {
      input.mu_prior = read_normal_prior(priors["mu"]);
    }
    read_vec_if_present(initial, "mu", input.initial.mu);
  }
}

template <typename Draws>
void read_draw_extensions(const Rcpp::List &object, Draws &draws)
{
  if (!has(object, "posterior")) {
    return;
  }
  const Rcpp::List posterior = object["posterior"];
  if constexpr (has_completion_draws<Draws>::value) {
    if (has(posterior, "y")) {
      read_draws_if_present(posterior["y"], "coeffs", draws.y);
    }
    if (has(posterior, "constraints_inv")) {
      read_draws_if_present(posterior["constraints_inv"], "coeffs", draws.constraints_inv);
    }
  }
  if constexpr (has_shrinkage_draws<Draws>::value) {
    if (has(posterior, "a") && Rf_isNewList(posterior["a"])) {
      const Rcpp::List a = posterior["a"];
      read_draws_if_present(a, "shrinkage", draws.a_shrinkage);
      read_draws_if_present(a, "local", draws.a_local);
    }
  }
  if constexpr (has_mu_draws<Draws>::value) {
    if (has(posterior, "mu")) {
      read_draws_if_present(posterior["mu"], "coeffs", draws.mu);
    }
  }
}

template <typename Draws>
void write_draw_extensions(Rcpp::List &posteriors, const Draws &draws)
{
  if constexpr (has_completion_draws<Draws>::value) {
    if (draws.y.n_elem > 0) {
      posteriors.push_back(Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.y)), "y");
    }
    if (draws.constraints_inv.n_elem > 0) {
      posteriors.push_back(
          Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.constraints_inv)),
          "constraints_inv");
    }
  }
  if constexpr (has_shrinkage_draws<Draws>::value) {
    if (draws.a_shrinkage.n_elem > 0 && posteriors.containsElementNamed("a") &&
        Rf_isNewList(posteriors["a"])) {
      Rcpp::List a(Rf_shallow_duplicate(posteriors["a"]));
      a.push_back(draws_to_r(draws.a_shrinkage), "shrinkage");
      if (draws.a_local.n_elem > 0) {
        a.push_back(draws_to_r(draws.a_local), "local");
      }
      posteriors["a"] = a;
    }
  }
  if constexpr (has_mu_draws<Draws>::value) {
    if (draws.mu.n_elem > 0) {
      posteriors.push_back(Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(draws.mu)), "mu");
    }
  }
}

} // namespace bayests_r

#endif // BVARTOOLS_BAYESTS_R_IO_H
