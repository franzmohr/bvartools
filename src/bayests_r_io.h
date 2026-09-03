#ifndef BVARTOOLS_BAYESTS_R_IO_H
#define BVARTOOLS_BAYESTS_R_IO_H

#include <RcppArmadillo.h>

#include "bayests/priors.h"
#include "bayests/spec.h"

#include <string>

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
  spec.p = optional_int(model, "p", 0);
  spec.m = optional_int(model, "m", 0);
  spec.s = optional_int(model, "s", 0);
  spec.n = optional_int(model, "n", 0);
  spec.n_restricted = optional_int(model, "n_restricted", 0);
  spec.rank = optional_int(model, "rank", 0);
  spec.k_beta = optional_int(model, "k_beta", 0);
  // Unobserved factors. Zero for every model this package has -- a dynamic
  // factor model is dfmtools' -- and read anyway, because the vendored core is
  // a whole mirror of upstream and its VarSpec carries the field.
  spec.n_factors = optional_int(model, "n_factors", 0);
  // Absent until add_forecast_input() has been called.
  spec.h = optional_int(model, "h", 0);
  spec.varsel = bayests::var_selection_from_string(optional_string(model, "varsel", "none"));
  spec.structural = optional_bool(model, "structural", false);

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

inline bayests::ConstantCointSpacePrior read_coint_space_prior_constant(const Rcpp::List &group)
{
  bayests::ConstantCointSpacePrior prior;
  read_double_if_present(group, "v_inv", prior.v_inv);
  read_mat_if_present(group, "p_tau_inv", prior.p_tau_inv);
  return prior;
}

/// The cointegration space prior of a model whose cointegration vectors move.
/// `rho` is the autoregression of their state equation and is not drawn, so it
/// arrives as a prior rather than as an initial value; a group that omits it
/// keeps the struct's default.
inline bayests::TvpCointSpacePrior read_coint_space_prior_tvp(const Rcpp::List &group)
{
  bayests::TvpCointSpacePrior prior;
  prior.initial_state = read_normal_prior(group);
  read_double_if_present(group, "rho", prior.rho);
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

} // namespace bayests_r

#endif // BVARTOOLS_BAYESTS_R_IO_H
